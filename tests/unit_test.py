import unittest
from signature_snvs import signature_snvs


class TestAnalysis(unittest.TestCase):
    def test_script(self):
        signature_snvs.signature_snvs_per_species.signature_snvs_per_species(
            species="Bacteroides_uniformis_57318",
            min_reads=5,
            start_index=1,
            end_index=200,
            config_file_path="/Users/leahbriscoe/Documents/FEASTX/Signature_SNVs/configs/config.yaml",
        )

    if __name__ == "__main__":
        unittest.main()

