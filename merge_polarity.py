import argparse
import networkx as nx
import requests
import os
import uuid
import pandas as pd

def merge(positive_graphml, negative_graphml, outputgraphml, RT_TOLERANCE=10, PPM_ERROR_TOLERANCE=20, masscolumn="precursor mass", rtcolumn="RTMean", output_summary_table=None, output_summary=None):
    #Loading Networks
    positive_G = nx.read_graphml(positive_graphml)
    negative_G = nx.read_graphml(negative_graphml)

    #Renaming
    positive_nodes_map = positive_G.nodes.data()
    negative_nodes_map = negative_G.nodes.data()

    renaming_positive_map = {}
    for positive_node in positive_nodes_map:
        positive_G.nodes[positive_node[0]]['ionization_mode'] = "positive"
        renaming_positive_map[positive_node[0]] = "positive-" + positive_node[0]

    renaming_negative_map = {}
    for negative_node in negative_nodes_map:
        negative_G.nodes[negative_node[0]]['ionization_mode'] = "negative"
        renaming_negative_map[negative_node[0]] = "negative-" + negative_node[0]
        
    positive_G = nx.relabel_nodes(positive_G, renaming_positive_map)
    negative_G = nx.relabel_nodes(negative_G, renaming_negative_map)


    positive_nodes_map = positive_G.nodes.data()
    negative_nodes_map = negative_G.nodes.data()

    plot_values = []

    RT_TOLERANCE = float(RT_TOLERANCE)
    PPM_ERROR_TOLERANCE = float(PPM_ERROR_TOLERANCE)

    new_edges = []
    for positive_node in positive_nodes_map:
        positive_id = positive_node[0]
        positive_mz = float(positive_node[1][masscolumn])
        positive_rt = float(positive_node[1][rtcolumn])
        
        for negative_node in negative_nodes_map:
            negative_id = negative_node[0]
            negative_mz = float(negative_node[1][masscolumn])
            negative_rt = float(negative_node[1][rtcolumn])
            
            mz_delta = positive_mz - negative_mz
            rt_delta = abs(positive_rt - negative_rt)
            
            expected_mz_delta = 1.007276 * 2
            
            mz_delta_error = abs(mz_delta - expected_mz_delta)
            mz_delta_error_ppm = mz_delta_error/negative_mz * 1000000
            
            if mz_delta < 1.9 or mz_delta > 2.1 or rt_delta > RT_TOLERANCE or mz_delta_error_ppm > PPM_ERROR_TOLERANCE:
                continue
                
            plot_values.append(rt_delta)
            new_edges.append([positive_id, negative_id])
            
    #Merging Networks
    merged_network = nx.compose(positive_G, negative_G)

    #Adding new edges
    node_dict = {}
    for node in list(merged_network.nodes(data=True)):
        node_dict[node[0]] = node[1]

    for edge in new_edges:
        node1 = node_dict[edge[0]]
        node2 = node_dict[edge[1]]

        mzdelta = abs(node1["precursor mass"] - node2["precursor mass"])
        rtdelta = abs(node1["RTMean"] - node2["RTMean"])

        try:
            smiles1 = node1["Smiles"]
            smiles2 = node2["Smiles"]
            
            if len(smiles1) > 5 and len(smiles2) > 5:
                similarity_url = "https://gnps-structure.ucsd.edu/structuresimilarity"
                r = requests.get(similarity_url, data={"smiles1": smiles1, "smiles2": smiles2})
                tanimoto = float(r.text)
            else:
                tanimoto = -1
        except:
            tanimoto = -1

        merged_network.add_edge(edge[0], edge[1], EdgeType="IonMode", mass_difference=mzdelta, rtdelta=rtdelta, tanimoto=tanimoto)

    #Evaluating the Network Merge
    all_edges = merged_network.edges(data=True)
    

    summary_list = []
    for edge in all_edges:
        edge_data = edge[2]
        node1 = node_dict[edge[0]]
        node2 = node_dict[edge[1]]
        
        if edge_data["EdgeType"] == "IonMode":
            try:
                output_dict = {}
                output_dict["scan1"] = edge[0]
                output_dict["scan2"] = edge[1]
                output_dict["smiles1"] = smiles1
                output_dict["smiles2"] = smiles2
                output_dict["mz1"] = node1["precursor mass"]
                output_dict["mz2"] = node2["precursor mass"]
                output_dict["mzdelta"] = edge_data["mass_difference"]
                output_dict["tanimoto"] = edge_data["tanimoto"]
                output_dict["rt1"] = node1["RTMean"]
                output_dict["rt2"] = node2["RTMean"]
                output_dict["rtdelta"] = edge_data["rtdelta"]

                summary_list.append(output_dict)
            except KeyboardInterrupt:
                raise
            except:
                continue

    #Saving out data
    nx.write_graphml(merged_network, outputgraphml)

    if output_summary_table is not None:
        pd.DataFrame(summary_list).to_csv(output_summary_table, sep="\t", index=False)

    total_pairs = len(summary_list)
    total_concordant_pairs = len([element for element in summary_list if element["tanimoto"] > 0.9 ])

    print("total_concordant_pairs", total_concordant_pairs)
    print("total_pairs", total_pairs)

    if output_summary is not None:
        with open(output_summary, "w") as output_summary_file:
            output_summary_file.write(f"total_concordant_pairs {total_concordant_pairs}\n")
            output_summary_file.write(f"total_pairs {total_pairs}\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--positive-network-task", dest="positivetask", default=None)
    parser.add_argument("--negative-network-task", dest="negativetask", default=None)
    parser.add_argument("--positive-graphml", dest="positivegraphml", default=None)
    parser.add_argument("--negative-graphml", dest="negativegraphml", default=None)
    parser.add_argument("--masscolumn", dest="masscolumn", default="precursor mass")
    parser.add_argument("--rtcolumn", dest="rtcolumn", default="RTMean")
    parser.add_argument("--output-graphml", dest="outputgraphml", default="merged_network.graphml")
    parser.add_argument("--output-summary-table", dest="outputsummarytable", default=None)
    parser.add_argument("--output-summary", dest="outputsummary", default=None)
    parser.add_argument("--RT_TOLERANCE", dest="RT_TOLERANCE", default="10")
    parser.add_argument("--PPM_ERROR_TOLERANCE", dest="PPM_ERROR_TOLERANCE", default="20")

    args = parser.parse_args()

    if args.positivegraphml is not None:
        positive_graphml = args.positivegraphml
    elif args.positivetask is not None:
        positive_url = "https://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task=%s&block=main&file=gnps_molecular_network_graphml/" % (args.positivetask)
        positive_graphml = str(uuid.uuid4()) + ".graphml"
        local_file = open(positive_graphml, "w")
        local_file.write(requests.get(positive_url).text)
        local_file.close()
    else:
        print("No input tasks")
        exit(1)

    if args.negativegraphml is not None:
        negative_graphml = args.negativegraphml
    elif args.negativetask is not None:
        negative_url = "https://gnps.ucsd.edu/ProteoSAFe/DownloadResultFile?task=%s&block=main&file=gnps_molecular_network_graphml/" % (args.negativetask)
        negative_graphml = str(uuid.uuid4()) + ".graphml"
        local_file = open(negative_graphml, "w")
        local_file.write(requests.get(negative_url).text)
        local_file.close()
    else:
        print("No input tasks")
        exit(1)

    merge(positive_graphml, \
        negative_graphml, \
        args.outputgraphml, RT_TOLERANCE=args.RT_TOLERANCE, \
        PPM_ERROR_TOLERANCE=args.PPM_ERROR_TOLERANCE, masscolumn=args.masscolumn, \
        rtcolumn=args.rtcolumn, output_summary_table=args.outputsummarytable, output_summary=args.outputsummary)




if __name__ == "__main__":
    main()
