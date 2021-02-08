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
        print(f"CalledProcessError exception occurred! Exception:\n{cpe_exp}")
        return cpe_exp
    except Exception as exp:
        print(f"Unknown exception occurred! Exception:\n{exp}")
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
    print(
        "Simple Python wrapper for MotEvo. "
        "For more details on usage check documentation at "
        "https://github.com/brlauuu/motevowrapper."
    )


def check_installation():
    result = shell_call(["motevo"])
    if result.returncode != 0:
        print(
            "MotEvo cannot be found on the system. Please follow the instructions "
            "at https://github.com/brlauuu/motevowrapper on how to download MotEvo."
        )
        return
    else:
        print(
            f"MotEvo successfully found on the system at: "
            f"{shell_call(['which', 'motevo'], verbose=True).stdout.decode('utf-8')}"
        )

    result = shell_call(["runUFE"])
    if result.returncode != 0:
        print(
            "runUFE executable cannot be found on the system. "
            "Please follow the instructions "
            "at https://github.com/brlauuu/motevowrapper on how to download MotEvo."
        )
    else:
        print(
            f"runUFE successfully found on the system at: "
            f"{shell_call(['which', 'runUFE'], verbose=True).stdout.decode('utf-8')}"
        )


def run_motevo(
    sequences_file=None,  # Or alignments file
    wm_path=None,
    working_directory="./",
    mode="TFBS",
    tree=None,
    ref_species=None,
    em_prior=None,
    ufe_wm_prior=None,
    ufe_wm_file=None,
    ufe_wm_len=None,
    background_prior=None,
    bgA=0.25,
    bgT=0.25,
    bgG=0.25,
    bgC=0.25,
    sites_file=None,
    priors_file=None,
    print_site_als=1,
    minposterior=0.1,
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

    if not sites_file:
        sites_file = f"sites_{pwm_name}"
    if not priors_file:
        priors_file = f"priors_{pwm_name}"

    # Load PWM length
    with open(wm_path, "r") as f:
        pwm_length = 0
        for line in f:
            if re.match("^\d+", line):
                pwm_length += 1

    # Create parameter file
    motevo_parameters_path = "motevo_parameters"
    with open(motevo_parameters_path, "w") as f:
        f.write(f"Mode {mode}\n")
        f.write(f"TREE {tree}\n")
        f.write(f"refspecies {ref_species}\n")
        f.write(f"EMprior {em_prior}\n")
        if ufe_wm_prior:
            f.write(f"UFEwmprior {ufe_wm_prior}\n")
        if ufe_wm_file:
            f.write(f"UFEwmfile {ufe_wm_file}\n")
        if ufe_wm_len and ufe_wm_len != "auto":
            f.write(f"UFEwmlen {ufe_wm_len}\n")
        elif ufe_wm_len == "auto":
            f.write(f"UFEwmlen {pwm_length}\n")
        f.write(f"bgprior {background_prior}\n")
        f.write(f"bg A {bgA}\n")
        f.write(f"bg T {bgT}\n")
        f.write(f"bg G {bgG}\n")
        f.write(f"bg C {bgC}\n")
        f.write(f"sitefile {sites_file}\n")
        f.write(f"priorfile {priors_file}\n")
        f.write(f"printsiteals {print_site_als}\n")
        f.write(f"minposterior {minposterior}\n")

    if verbose:
        print(f"Generated parameters file at: {motevo_parameters_path}")

    # Remove existing MotEvo outputs
    if os.path.exists(sites_file):
        os.remove(sites_file)
    if os.path.exists(priors_file):
        os.remove(priors_file)

    # Run MotEvo
    result = shell_call(["motevo", sequences_file, motevo_parameters_path, wm_path])

    # Check result
    if result.returncode == 0:
        if verbose:
            print(
                f"MotEvo ran successfully! Please"
                f"check results at: {sites_file} and {priors_file}"
            )
    else:
        print("MotEvo run failed!")

    # Change back to working directory
    os.chdir(cwd)

    return (
        os.path.join(working_directory, sites_file),
        os.path.join(working_directory, priors_file),
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
