import pathlib
import unittest

import pandas as pd

from map_tree import make_abs_pres_networkx

class TestAbsencePresence(unittest.TestCase):
    def setUp(self) -> None:
        self.base_dir = "/home/kat/Documents/work/DTU/bioengineering10/test_data/"

    def test_map_families(self):
        clustering_file = pathlib.Path(self.base_dir) / "success_test" / "Network_Annotations_Full_success_test_mibig"
        clusters_to_find = {"Test_organism_1": {"testcluster1", "testcluster3", "testcluster5", "testcluster7",
                                                "testcluster8"},
                            "Test_organism_2": {"testcluster2", "testcluster4", "testcluster6"}}
        bgcs_to_find = {"BGC_1": "Testomycin"}
        types_to_find = {"testcluster1": "NRPS", "testcluster2": "PKS", 
                         "testcluster3": "other", "testcluster4": "NRPS",
                         "testcluster5": "other", "testcluster6": "other",
                         "testcluster7": "other", "testcluster8": "PKS",
                         "BGC_1": "other"}
        clusters_found, bgcs_found, types_found = make_abs_pres_networkx.map_clusters_to_names(clustering_file)

        assert clusters_found == clusters_to_find
        assert bgcs_found == bgcs_to_find
        assert types_found == types_to_find

    def test_connected_from_df(self):
        network_file = pathlib.Path(self.base_dir) / "success_test" / "full_network.tsv"
        test_network = pd.read_csv(network_file, sep="\t")
        test_types = {"Cluster_1": "test", "Cluster_2": "test", "Cluster_3": "test", "Cluster_4": "test",
                      "Cluster_5": "test"}
        test_components = make_abs_pres_networkx.get_connected_from_df(test_network, test_types)
        assert len(test_components) == 2
        assert test_components[0].compound_class == "test"
        assert test_components[1].compound_class == "test"

    def test_get_families(self):
        network_file = pathlib.Path(self.base_dir) / "success_test" / "full_network.tsv"
        test_network = pd.read_csv(network_file, sep="\t")
        test_types = {"Cluster_1": "test", "Cluster_2": "test", "Cluster_3": "test", "Cluster_4": "test",
                      "Cluster_5": "test"}
        test_components = make_abs_pres_networkx.get_connected_from_df(test_network, test_types)
        bgc_fakes = {"BGC_1": "Testomcyin"}
        test_families = make_abs_pres_networkx.get_families(test_components, bgc_fakes)
        assert len(test_families["1_test"]) == 3
        assert len(test_families["2_test"]) == 2
        assert "Cluster_1" in test_families["1_test"]
        assert "Cluster_5" not in test_families["1_test"]
        assert "Cluster_5" in test_families["2_test"]

    def test_merge_networks(self):
        result_dirs = []
        success_dir = pathlib.Path(self.base_dir) / "success_test"
        for entry in success_dir.iterdir():
            if entry.is_dir():
                result_dirs.append(entry)
        test_networks = make_abs_pres_networkx.merge_all_bigscape_networks(result_dirs)
        # do we have everything?
        assert test_networks[["Clustername 1"]].size == 3
        # is it in the right place?
        assert "Cluster_1" in test_networks[["Clustername 1"]].values
        # did (roughly) averaging the duplicates work?
        test_mean_for_duplicates = test_networks[test_networks['Clustername 1'] == 'Cluster_1'][['Raw distance']]
        assert 0.35 in test_mean_for_duplicates.values

    def test_matrix(self):
        fake_families = {"1_test": {"Cluster_1", "Cluster_2", "Cluster_3"}, "2_test": {"Cluster_4", "Cluster_5"}}
        fake_genomes = {"Test_org_1": {"Cluster_1"}, "Test_org_2": {"Cluster_3", "Cluster_4"},
                        "Test_org_3": {"Cluster_5"}}
        true_matrix = "#NAMES,1_test,2_test\nTest_org_1,1,0\nTest_org_2,1,1\nTest_org_3,0,1\n"
        test_matrix = make_abs_pres_networkx.get_absence_presence_matrix(fake_families, fake_genomes)
        assert test_matrix == true_matrix

    def test_preordered_matrix(self):
        fake_families = {"1_test": {"Cluster_1", "Cluster_2", "Cluster_3"}, "2_test": {"Cluster_4", "Cluster_5"}}
        fake_genomes = {"Test_org_1": {"Cluster_1"}, "Test_org_2": {"Cluster_3", "Cluster_4"},
                        "Test_org_3": {"Cluster_5"}}
        fake_order = ["2", "1"]
        true_matrix = "#NAMES,2_test,1_test\nTest_org_1,0,1\nTest_org_2,1,1\nTest_org_3,1,0\n"
        test_matrix = make_abs_pres_networkx.get_preordered_matrix(fake_families, fake_genomes, fake_order)
        assert test_matrix == true_matrix

    def test_extract_order(self):
        order_file = pathlib.Path(self.base_dir) / "test_order"
        true_order = ["5", "4", "1", "2", "3"]
        extract_order = make_abs_pres_networkx.parse_family_order(order_file)
        assert extract_order == true_order