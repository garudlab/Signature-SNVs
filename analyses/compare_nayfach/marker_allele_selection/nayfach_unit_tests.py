import unittest

from nayfach_helpers import get_column_names, read_two_files, check_marker_allele


class TestAnalysis(unittest.TestCase):
    # def setUp(self):

    # ref allele is the marker allele in C
    def check_marker_alleles_1(self):
        marker_alleles = check_marker_allele(
            focal_sample_names=["A", "B", "C"],
            focal_depth=[10, 50, 60],
            focal_freq=[0.01, 0.01, 0.2],
            background_samples=["D", "E", "F"],
            background_depth=[30, 30, 30],
            background_freq=[0.01, 0.01, 0.05],
        )
        print(marker_alleles)
        self.assertEqual(marker_alleles, ["C"])

    #  ref is the marker allele in all

    def test_check_marker_alleles_2(self):
        marker_alleles = check_marker_allele(
            focal_sample_names=["A", "B", "C"],
            focal_depth=[33, 50, 60],
            focal_freq=[0.11, 0.10, 0.90],
            background_samples=["D", "E", "F"],
            background_depth=[30, 30, 30],
            background_freq=[0.01, 0.01, 0.05],
        )
        print(marker_alleles)
        self.assertEqual(sorted(marker_alleles), ["A", "B", "C"])

    # marker alleles in different ways. First two are from alt, second two are from ref
    def test_check_marker_alleles_3(self):
        marker_alleles = check_marker_allele(
            focal_sample_names=["A", "B", "C", "D"],
            focal_depth=[33, 50, 60, 40],
            focal_freq=[0.05, 0.01, 0.90, 0.90],
            background_samples=["G", "E", "F"],
            background_depth=[20, 20, 10],
            background_freq=[0.01, 0.01, 0.05],
        )
        print(marker_alleles)
        self.assertEqual(sorted(marker_alleles), ["A", "B", "C", "D"])


if __name__ == "__main__":
    unittest.main()
