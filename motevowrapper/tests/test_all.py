import unittest
import os
import pandas as pd
from pandas.util.testing import assert_equal, assert_frame_equal

from motevowrapper.motevowrapper import (
    parse_sites,
    parse_priors,
    run_motevo,
    run_ufe,
    shell_call,
)

BASE_PATH = os.path.abspath(os.path.dirname(__file__))
DATA_PATH = os.path.join(BASE_PATH, "data")
OUTPUT_PATH = os.path.join(BASE_PATH, "output")
if not os.path.exists(OUTPUT_PATH):
    os.mkdir(OUTPUT_PATH)


class TestMotevoWrapper(unittest.TestCase):
    def test_parsing_sites(self):
        motifs = ["REST"] * 15
        reference_promoters = [
            "danRer11_chr21_11468142_11469142_+",
            "danRer11_chr10_21587502_21588502_+",
            "danRer11_chr6_43091675_43092675_-",
            "danRer11_chr5_29749877_29750877_-",
            "danRer11_chr4_17390591_17391591_-",
            "danRer11_chr21_6113805_6114805_+",
            "danRer11_chr25_5034843_5035843_+",
            "danRer11_chr25_5034843_5035843_+",
            "danRer11_chr9_68434_69447_-",
            "danRer11_chr8_51404306_51405306_-",
            "danRer11_chr8_51404306_51405306_-",
            "danRer11_chr8_51404306_51405306_-",
            "danRer11_chr8_51404306_51405306_-",
            "danRer11_chr8_51404306_51405306_-",
            "danRer11_chr17_12384808_12385808_-",
        ]

        reference_binding_strands = [
            "-",
            "-",
            "-",
            "+",
            "-",
            "+",
            "+",
            "+",
            "+",
            "-",
            "-",
            "-",
            "-",
            "-",
            "-",
        ]

        motif_coordinates = [
            "471-491",
            "75-95",
            "869-889",
            "602-622",
            "182-202",
            "936-956",
            "896-916",
            "896-916",
            "671-691",
            "742-762",
            "742-762",
            "742-762",
            "742-762",
            "742-762",
            "776-796",
        ]

        posteriors = [
            0.546355,
            0.991828,
            0.362056,
            0.999812,
            0.997489,
            0.107395,
            0.999404,
            0.999404,
            0.319236,
            0.716461,
            0.716461,
            0.716461,
            0.716461,
            0.716461,
            0.999844,
        ]

        aligned_promoters = [
            "danRer11_chr21_11468142_11469142_+",
            "danRer11_chr10_21587502_21588502_+",
            "danRer11_chr6_43091675_43092675_-",
            "danRer11_chr5_29749877_29750877_-",
            "danRer11_chr4_17390591_17391591_-",
            "danRer11_chr21_6113805_6114805_+",
            "danRer11_chr25_5034843_5035843_+",
            "ictPun_chr14_23588333_23589324_-",
            "danRer11_chr9_68434_69447_-",
            "danRer11_chr8_51404306_51405306_-",
            "astMex_chr6_11399458_11400473_-",
            "pygNat_KV575330.1_1811789_1812923_-",
            "ictPun_chr5_4628508_4629508_-",
            "esoLuc_chrLG13_14038768_14039768_+",
            "danRer11_chr17_12384808_12385808_-",
        ]

        scores = [
            20.3009,
            13.495,
            9.07795,
            18.1038,
            16.0274,
            8.6284,
            18.1535,
            14.0058,
            9.37126,
            5.93671,
            16.4049,
            16.4049,
            12.4296,
            10.2324,
            18.9164,
        ]

        binding_sequences = [
            "AGCGCTGTCCTTGGTGCTGAC",
            "CTCGTTGTCCAAGGTGCTGAA",
            "TTATTTGTCCATGGTTCTGAT",
            "GTACCTGTCCTTGGTGCTGAA",
            "GCTGCTCTCCAAGGTACTGAA",
            "GCTGCTGTTCCGCGTGCTGGA",
            "ATCGCTGTCCATGGTGCTGCA",
            "GCTACTGTCCATGGTGCTGTT",
            "CCCGCTGTCCGCCGTTCTGGA",
            "GGCGCTGTCTTTAGTACAGGA",
            "GGCGCTGTCCTTGGTGCAGGA",
            "GGCGCTGTCCTTGGTGCAGGA",
            "GGCGCTGTCTTTGGTGCAGGA",
            "GGCACTGTCTTTGGTGCAGGA",
            "AGCGCTCTCCGCGGTGCTGAA",
        ]

        dict = {
            "motif": motifs,
            "reference_promoter": reference_promoters,
            "reference_binding_strand": reference_binding_strands,
            "motif_coordinates": motif_coordinates,
            "posterior": posteriors,
            "aligned_promoter": aligned_promoters,
            "score": scores,
            "binding_sequence": binding_sequences,
        }

        values_df = pd.DataFrame.from_dict(dict)
        results_df = parse_sites(os.path.join(DATA_PATH, "sites_REST.wm"))

        assert_frame_equal(values_df, results_df, check_dtype=False)

    def test_parsing_priors(self):
        motifs = ["REST", "background", "UFEwm"]
        final_priors = [0.00310981, 0.828626, 0.168265]
        nr_of_sites = [7.04002, 1875.85, 380.919]
        densities = [0.0147501, 0.187155, 0.798095]

        dict = {
            "motif": motifs,
            "final_prior": final_priors,
            "nr_of_sites": nr_of_sites,
            "density": densities,
        }

        values_df = pd.DataFrame.from_dict(dict)
        results_df = parse_priors(os.path.join(DATA_PATH, "priors_REST.wm"))

        assert_frame_equal(values_df, results_df, check_dtype=False)

    def test_motevo_run(self):
        result = run_motevo(
            sequences_file=os.path.join(DATA_PATH, "zebrafish_alignments.aln"),
            working_directory=OUTPUT_PATH,
            wm_path=os.path.join(DATA_PATH, "pwmdir", "REST.wm"),
            tree="((((astMex:0.415917,pygNat:0.449133):0.099801,ictPun:0.50305):0.04815395,danRer11:0.55291):0.0098669,esoLuc:0.7121605);",
            ref_species="danRer11",
            em_prior=0,
            ufe_wm_prior=500,
            ufe_wm_file=os.path.join(DATA_PATH, "UFEmodel"),
            ufe_wm_len="auto",
            background_prior=0.8,
        )
        self.assertEqual(result[0], os.path.join(OUTPUT_PATH, "sites_REST.wm"))
        results_1 = parse_sites(os.path.join(DATA_PATH, "sites_REST.wm"))
        results_2 = parse_sites(os.path.join(OUTPUT_PATH, "sites_REST.wm"))
        assert_frame_equal(results_1, results_2, check_dtype=False)

        self.assertEqual(result[1], os.path.join(OUTPUT_PATH, "priors_REST.wm"))
        results_1 = parse_priors(os.path.join(DATA_PATH, "priors_REST.wm"))
        results_2 = parse_priors(os.path.join(OUTPUT_PATH, "priors_REST.wm"))
        assert_frame_equal(results_1, results_2, check_dtype=False)

    def test_installation(self):
        result = shell_call(["motevo"])
        self.assertEqual(result.returncode, 0)
        result = shell_call(["runUFE"])
        self.assertEqual(result.returncode, 1)

    def test_ufe_run(self):
        run_ufe(
            tree_file_path=os.path.join(DATA_PATH, "tree_file"),
            output_path=os.path.join(OUTPUT_PATH, "UFEmodel"),
        )

        f = open(os.path.join(OUTPUT_PATH, "UFEmodel"))
        result = f.readlines()
        f.close()

        g = open(os.path.join(DATA_PATH, "UFEmodel"))
        check = g.readlines()
        g.close()

        self.assertEqual(len(result), len(check))

        for i in range(len(result)):
            self.assertEqual(result[i], check[i])


if __name__ == "__main__":
    unittest.main()
