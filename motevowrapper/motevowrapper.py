import os
import re
import logging
import subprocess
import pandas as pd


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
        err = cpe_exp.stderr.decode("utf-8")
        print(f"CalledProcessError exception occurred! Error: {err}")
        print(f"Full exception:\n{cpe_exp}")
        return cpe_exp
    except Exception as exp:
        err = exp.stderr.decode("utf-8")
        print(f"Unknown exception occurred! Error:{err}")
        print(f"Full exception:\n{exp}")
        return exp


def parse_sites(path, verbose=False):
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
                # if j == len(lines):
                #     break
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


def parse_priors(path, verbose=False):
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
        "Simple Python wrapper for MotEvo."
        "For more details, go to https://github.com/brlauuu/motevowrapper"
    )


def check_installation():
    result = shell_call(["motevo"])
    if result.returncode != 0:
        print(
            "MotEvo cannot be found on the system. Please follow the instructions"
            "at https://github.com/brlauuu/motevowrapper on how to download MotEvo."
        )
    else:
        print(
            f"MotEvo successfully found on the system at: "
            f"{shell_call(['which', 'motevo'], verbose=True).stdout.decode('utf-8')}"
        )


def run():
    assert False, "Not implemented yet!"
