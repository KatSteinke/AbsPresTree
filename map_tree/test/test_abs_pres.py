import pathlib
import unittest

from map_tree import make_abs_pres


class TestAbsencePresence(unittest.TestCase):
    def setUp(self) -> None:
        self.base_dir = "/home/katismash/Documents/work/DTU/bioengineering10/test_data/"

    def test_get_families(self):
        families_file = pathlib.Path(self.base_dir) / "success_test" / "NRPS" / "NRPS_clustering_success.tsv"
        families = {}
        families_to_find = {"1": {"testcluster1", "testcluster4"}}
        make_abs_pres.get_families(families_file, families)

        assert families == families_to_find

    def test_match_clusters(self):
        clustering_file = pathlib.Path(self.base_dir) / "success_test" / "Network_Annotations_Full_success_test"
        clusters_to_find = {"Test_organism_1": {"testcluster1", "testcluster3", "testcluster5", "testcluster7",
                                                "testcluster8"},
                            "Test_organism_2": {"testcluster2", "testcluster4", "testcluster6"}}
        clusters_found = make_abs_pres.match_clusters_to_genomes(clustering_file)

        assert clusters_found == clusters_to_find

    def test_make_matrix(self):
        success_dir = self.base_dir + "success_test/"
        success_matrix = "#NAMES\t1\t2\t3\t4\t5\nTest_organism_1\t2\t0\t1\t1\t1\nTest_organism_2\t1\t1\t0\t1\t0"
        result_matrix = make_abs_pres.bigscape_to_matrix(success_dir)
        assert success_matrix == result_matrix

    def test_catch_missing_files(self):
        no_summary = self.base_dir + "fail_test_no_network"
        with self.assertRaisesRegex(IOError, "No summary file found."):
            make_abs_pres.bigscape_to_matrix(no_summary)

        no_dirs = self.base_dir + "fail_test_no_folders"
        with self.assertRaisesRegex(IOError, "No results found."):
            make_abs_pres.bigscape_to_matrix(no_dirs)

        no_clusters = self.base_dir + "fail_test_no_clusters"
        with self.assertRaisesRegex(ValueError, "No clustering files in directory other"):
            make_abs_pres.bigscape_to_matrix(no_clusters)

    def test_duplicate_file(self):
        duplicate_summary = self.base_dir + "fail_test_duplicates"
        with self.assertRaisesRegex(ValueError, "Multiple summary files found."):
            make_abs_pres.bigscape_to_matrix(duplicate_summary)

        duplicate_clusters = self.base_dir + "fail_test_multiple_clusters"
        with self.assertRaisesRegex(ValueError, "Multiple clustering files in directory other"):
            make_abs_pres.bigscape_to_matrix(duplicate_clusters)

    def test_catch_summary_format_fails(self):
        summary_no_header = pathlib.Path(self.base_dir) / "fail_test" / "Network_Annotations_Full_fail_no_header"
        with self.assertRaisesRegex(ValueError, "Header line of summary file missing."):
            make_abs_pres.match_clusters_to_genomes(summary_no_header)

    def catch_cluster_format_fails(self):
        families = {}

        cluster_no_header = pathlib.Path(self.base_dir) / "fail_test" / "NRPS" / "NRPS_clustering_no_header.tsv"
        with self.assertRaisesRegex(ValueError, "Header line of file NRPS_clustering_no_header.tsv missing"):
            make_abs_pres.get_families(cluster_no_header, families)

        cluster_wrong_columns = pathlib.Path(self.base_dir) / "fail_test" / "other" / "other_clustering_column_fail.tsv"
        with self.assertRaisesRegex(ValueError, "Cluster file must contain two columns. File other_clustering_column_fail.tsv contains 1."):
            make_abs_pres.get_families(cluster_wrong_columns, families)
