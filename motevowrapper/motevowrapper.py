import os
import re
import logging
import subprocess
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


logger = logging.getLogger(__name__)


def shell_call(command, verbose=False):
    """
    Method that performs a shell call while not forwarding stdout and stderr.
    """
    try:
        if verbose:
            result = subprocess.run(command, capture_output=True)
        else:
            result = subprocess.run(
                command, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            )
        return result

    except subprocess.CalledProcessError as cpe_exp:
        logger.error(f"CalledProcessError exception occurred! Exception:\n{cpe_exp}")
        return cpe_exp
    except Exception as exp:
        logger.error(f"Unknown exception occurred! Exception:\n{exp}")
        return exp


def parse_sites(path):
    if not os.path.exists(path):
        logger.error(f"Path doesn't exist: {path}")

    prep = []
    with open(path, "r") as f:
        lines = f.readlines()
        i = 0
        while i < len(lines):
            if re.match(r"^\d+", lines[i]):
                columns = lines[i].strip().split(" ")
                motif_coordinates = columns[0].strip()
                reference_binding_strand = columns[1].strip()
                posterior = float(columns[2])
                motif = columns[3].strip()
                reference_promoter = columns[4].strip()

                j = i + 1
                while j < len(lines) and not re.match(r"^\d+", lines[j]):
                    columns = lines[j].strip().split(" ")
                    sequence = columns[0]
                    score = float(columns[1])
                    aligned_promoter = columns[2]

                    prep.append(
                        {
                            "motif": motif,
                            "reference_promoter": reference_promoter,
                            "reference_binding_strand": reference_binding_strand,
                            "motif_coordinates": motif_coordinates,
                            "posterior": posterior,
                            "aligned_promoter": aligned_promoter,
                            "score": score,
                            "binding_sequence": sequence,
                        }
                    )
                    j += 1

                i = j
            else:
                logger.error(f"Missing sequence line at line {i} in {path}!")
                break

    return pd.DataFrame(prep)


def parse_priors(path):
    if not os.path.exists(path):
        logger.error(f"Path doesn't exist: {path}")

    prep = []
    with open(path, "r") as f:
        f.readline()  # header

        try:
            motif_results = f.readline().strip().split(" ")
            prep.append(
                {
                    "motif": motif_results[0],
                    "final_prior": float(motif_results[1]),
                    "nr_of_sites": float(motif_results[2]),
                    "density": float(motif_results[3]),
                }
            )

            background_results = f.readline().strip().split(" ")
            prep.append(
                {
                    "motif": background_results[0],
                    "final_prior": float(background_results[1]),
                    "nr_of_sites": float(background_results[2]),
                    "density": float(background_results[3]),
                }
            )

            ufe_results = f.readline().strip().split(" ")
            prep.append(
                {
                    "motif": ufe_results[0],
                    "final_prior": float(ufe_results[1]),
                    "nr_of_sites": float(ufe_results[2]),
                    "density": float(ufe_results[3]),
                }
            )
        except:
            logger.error(f"Not enough lines in {path}")
            return None

    return pd.DataFrame(prep)


def print_help():
    logger.info(
        "Simple Python wrapper for MotEvo. "
        "For more details on usage check documentation at "
        "https://github.com/brlauuu/motevowrapper."
    )


def check_installation():
    result = shell_call(["motevo"])
    if result.returncode != 0:
        logger.error(
            "MotEvo cannot be found on the system. Please follow the instructions "
            "at https://github.com/brlauuu/motevowrapper on how to download MotEvo."
        )
        return
    else:
        logger.info(
            f"MotEvo successfully found on the system at: "
            f"{shell_call(['which', 'motevo'], verbose=True).stdout.decode('utf-8')}"
        )

    result = shell_call(["runUFE"])
    if result.returncode != 1:
        logger.error(
            "runUFE executable cannot be found on the system. "
            "Please follow the instructions "
            "at https://github.com/brlauuu/motevowrapper on how to download MotEvo."
        )
    else:
        logger.info(
            f"runUFE successfully found on the system at: "
            f"{shell_call(['which', 'runUFE'], verbose=True).stdout.decode('utf-8')}"
        )


def run_motevo(
    sequences_file=None,
    wm_path=None,
    working_directory="./",
    Mode="TFBS",
    refspecies=None,
    TREE=None,
    restrictparses=None,
    singlestrand=None,
    bgprior=None,
    EMprior=0,
    priordiff=None,
    UFEwmprior=None,
    markovorderBG=None,
    bgA=0.25,
    bgT=0.25,
    bgG=0.25,
    bgC=0.25,
    mybgfile=None,
    UFEwmfile=None,
    UFEwmlen=None,
    UFEprint=None,
    UFEwmproffile=None,
    sitefile=None,
    priorfile=None,
    loglikfile=None,
    minposterior=0.1,
    printsiteals=1,
    minposteriorWM=None,
    wmdiff=None,
    CRMfile=None,
    winlen=None,
    steplen=None,
    try_until_succeeding=False,
    verbose=False,
):
    # Check if MotEvo is installed
    assert shell_call(["motevo"]).returncode == 0, (
        "Could not find MotEvo. Please check installation"
        "first by running `check_installation()` method!"
    )

    # Change directory to working_directory
    cwd = os.getcwd()
    os.chdir(working_directory)

    # Read Position Weight Matrix (PWM) name
    pwm_name = wm_path[wm_path.rfind("/") + 1 :]

    if not sitefile:
        sitefile = f"sites_{pwm_name}"
    if not priorfile:
        priorfile = f"priors_{pwm_name}"

    if not TREE:
        TREE = f"({refspecies}: 1.0);"

    # Load PWM length
    with open(wm_path, "r") as f:
        pwm_length = 0
        for line in f:
            if re.match("^\d+", line):
                pwm_length += 1

    # Create parameter file
    motevo_parameters_path = "motevo_parameters"
    with open(motevo_parameters_path, "w") as f:
        if Mode:
            f.write(f"Mode {Mode}\n")
        if TREE:
            f.write(f"TREE {TREE}\n")
        if refspecies:
            f.write(f"refspecies {refspecies}\n")
        if bgprior:
            f.write(f"bgprior {bgprior}\n")
        if minposterior:
            f.write(f"minposterior {minposterior}\n")
        if sitefile:
            f.write(f"sitefile {sitefile}\n")
        if priorfile:
            f.write(f"priorfile {priorfile}\n")
        if bgA:
            f.write(f"bg A {bgA}\n")
        if bgT:
            f.write(f"bg T {bgT}\n")
        if bgG:
            f.write(f"bg G {bgG}\n")
        if bgC:
            f.write(f"bg C {bgC}\n")
        if mybgfile:
            f.write(f"mybgfile {mybgfile}\n")
        if EMprior:
            f.write(f"EMprior {EMprior}\n")
        if minposteriorWM:
            f.write(f"minposteriorWM {minposteriorWM}\n")
        if UFEwmprior:
            f.write(f"UFEwmprior {UFEwmprior}\n")
        if UFEwmfile:
            f.write(f"UFEwmfile {UFEwmfile}\n")
        if UFEwmlen and UFEwmlen != "auto":
            f.write(f"UFEwmlen {UFEwmlen}\n")
        elif UFEwmlen == "auto":
            f.write(f"UFEwmlen {pwm_length}\n")
        if UFEwmproffile:
            f.write(f"UFEwmproffile {UFEwmproffile}\n")
        if UFEprint:
            f.write(f"UFEprint {UFEprint}\n")
        if wmdiff:
            f.write(f"wmdiff {wmdiff}\n")
        if winlen:
            f.write(f"winlen {winlen}\n")
        if steplen:
            f.write(f"steplen {steplen}\n")
        if markovorderBG:
            f.write(f"markovorderBG {markovorderBG}\n")
        if priordiff:
            f.write(f"priordiff {priordiff}\n")
        if restrictparses:
            f.write(f"restrictparses {restrictparses}\n")
        if loglikfile:
            f.write(f"loglikfile {loglikfile}\n")
        if CRMfile:
            f.write(f"CRMfile {CRMfile}\n")
        if singlestrand:
            f.write(f"singlestrand {singlestrand}\n")
        if printsiteals:
            f.write(f"printsiteals {printsiteals}\n")

    if verbose:
        logger.info(f"Generated parameters file at: {motevo_parameters_path}.")

    # Remove existing MotEvo outputs
    if os.path.exists(sitefile):
        os.remove(sitefile)
    if os.path.exists(priorfile):
        os.remove(priorfile)

    # Setting the status of running MotEvo
    status = False

    while not status:
        command = ["motevo", sequences_file, motevo_parameters_path, wm_path]

        if verbose:
            print(f"MotEvo shell command:\n" f"{' '.join(command)}")

        # Run MotEvo
        result = shell_call(command, verbose=True)

        # Writing motevo report
        with open("motevo_report", "w") as f:
            f.write(result.stdout.decode("utf-8"))

        # Check result
        if result.returncode == 0:
            if verbose:
                logger.info(
                    f"MotEvo ran successfully! Please"
                    f"check results at: {sitefile} and {priorfile}.\n"
                    f"Check report at motevo_report."
                )
            status = True
        else:
            logger.error("MotEvo run failed!")
            status = False

        # Check if files were generated
        if not os.path.exists(sitefile):
            logger.error("MotEvo did not generate sites file.")
            status = False

        if not os.path.exists(priorfile):
            logger.error("MotEvo did not generate priors file.")
            status = False

        # In case user wants to run only once, we break the loop
        # Otherwise, MotEvo will be attempted to run until it succeeds
        if not try_until_succeeding:
            break

    # Change back to working directory
    os.chdir(cwd)

    return (
        os.path.join(working_directory, sitefile),
        os.path.join(working_directory, priorfile),
    )


def run_ufe(
    tree_file_path, bg_A=0.25, bg_C=0.25, bg_G=0.25, bg_T=0.25, output_path=None
):
    # Check if runUFE is installed
    assert shell_call(["runUFE"]).returncode == 1, (
        "Could not find runUFE. Please check installation"
        "first by running `check_installation()` method!"
    )

    if not output_path:
        output_path = "UFE_model"

    result = shell_call(
        ["runUFE", tree_file_path, str(bg_A), str(bg_C), str(bg_G), str(bg_T)],
        verbose=True,
    )

    with open(output_path, "w") as f:
        f.write(result.stdout.decode("utf-8"))

    return os.path.join(os.getcwd(), output_path)


def plot_site_distribution(motif, df, kind="ecdf"):
    sns.set_context("talk")
    sns.set_style("whitegrid")
    sns.displot(data=df.groupby("reference_promoter").sum(), kind=kind, x="posterior")
    plt.title(motif)
    plt.show()
