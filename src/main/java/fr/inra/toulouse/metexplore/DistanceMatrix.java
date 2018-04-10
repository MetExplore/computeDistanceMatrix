package fr.inra.toulouse.metexplore;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

import fr.inra.toulouse.metexplore.met4j_core.biodata.BioChemicalReaction;
import fr.inra.toulouse.metexplore.met4j_core.biodata.BioNetwork;
import fr.inra.toulouse.metexplore.met4j_core.biodata.BioPhysicalEntity;
import fr.inra.toulouse.metexplore.met4j_core.io.Sbml2BioNetworkLite;
import fr.inra.toulouse.metexplore.met4j_graph.computation.algo.ShortestPath;
import fr.inra.toulouse.metexplore.met4j_graph.computation.transform.GraphFilter;
import fr.inra.toulouse.metexplore.met4j_graph.computation.weighting.DefaultWeightPolicy;
import fr.inra.toulouse.metexplore.met4j_graph.computation.weighting.DegreeWeightPolicy;
import fr.inra.toulouse.metexplore.met4j_graph.computation.weighting.WeightUtils;
import fr.inra.toulouse.metexplore.met4j_graph.computation.weighting.WeightsFromFile;
import fr.inra.toulouse.metexplore.met4j_graph.core.BioPath;
import fr.inra.toulouse.metexplore.met4j_graph.core.WeightingPolicy;
import fr.inra.toulouse.metexplore.met4j_graph.core.compound.CompoundGraph;
import fr.inra.toulouse.metexplore.met4j_graph.core.compound.ReactionEdge;
import fr.inra.toulouse.metexplore.met4j_graph.io.Bionetwork2BioGraph;

public class DistanceMatrix {

	/**
	 * Example args:
	 * E:\Recherche\RESEARCH-DB-Networks\Homo_sapiens\recon2.v03_ext_noCompartment_noTransport.xml E:\workspace\NetworkComparison\src\main\resources\DistanceMatrix\Sonia\fingerprint.txt E:\workspace\NetworkComparison\src\main\resources\DistanceMatrix\Sonia\distanceMatrix.txt
	 * @param args
	 */
	
	public static void main(String[] args) {
		
		System.err.println("----GO");

		/**
		 * Path to SBML network, for example E:\Recherche\RESEARCH-DB-Networks\Homo sapiens\recon2.v03_ext_noCompartment_noTransport.xml
		 */
		String sbmlPath=args[0];
		/**
		 * Graph creation
		 */
		//import the network
		Sbml2BioNetworkLite importer=new Sbml2BioNetworkLite(sbmlPath,true);
		importer.setNotesValueSeparator(" || ");
		//Create the BioNetwork
		BioNetwork net=importer.getBioNetwork();
		//Turn the BioNetwork into a compound graph
		CompoundGraph compoundGraph=(new Bionetwork2BioGraph(net)).getCompoundGraph();
		
		/***
		 * Creating weight on edges. Can be done using lightest path or importing weight from a text file like atom mapping weights
		 */
//		//---------------Import weights on edges proportion of carbon transfer
//		WeightsFromFile<BioPhysicalEntity, ReactionEdge, CompoundGraph> weights=new WeightsFromFile<BioPhysicalEntity, ReactionEdge, CompoundGraph>(args[3],true);
//		weights.setWeight(compoundGraph);
//		GraphFilter.weightFilter(compoundGraph,0,GraphFilter.EQUALITY);
//		WeightUtils.removeEdgeWithNaNWeight(compoundGraph);
//		
//		DefaultWeightPolicy<BioPhysicalEntity, ReactionEdge, CompoundGraph> wpolicy=new DefaultWeightPolicy<BioPhysicalEntity, ReactionEdge, CompoundGraph>();
//		wpolicy.setWeight(compoundGraph);
		
		
		// TODO Auto-generated method stub
		/**
		 * Get metabolic profile as a tabular file with one column with metabolite identifiers used in the network
		 */
		HashSet<BioPhysicalEntity> metabolites=new HashSet<BioPhysicalEntity>();
		try {
			BufferedReader reader=new BufferedReader(new FileReader(args[1]));
			try {
				String line=reader.readLine();
				while(line!=null)
				{
					String[] elements=line.split("\t");
					BioPhysicalEntity metabolite=net.getBioPhysicalEntityById(elements[0]);
					metabolites.add(metabolite);
					if (metabolite==null)
						System.out.println(elements[0]+" not found ");
					line=reader.readLine();
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		} catch (FileNotFoundException e) {
			System.err.println("Please provide fingerprint file");
			e.printStackTrace();
		}
		

		/**
		 * Main part of the function that will create a distance matrix and print it out in a file.
		 */
		try {
			BufferedWriter writer=new BufferedWriter(new FileWriter(args[2]));
			
			/**
			 * Weight based on node degree squared. Will correspond to lightest path square
			 */
//			DegreeWeightPolicy weighting=new DegreeWeightPolicy(2);
//			weighting.setWeight(compoundGraph);
			
			/**
			 * Import weight from atom mapping set up in file
			 * TODO: add this file as an optional argument of the function and make it generic so any weighting can be imported
			 */
			compoundGraph.removeVertex(compoundGraph.getVertex("M_co2"));
			WeightsFromFile<BioPhysicalEntity, ReactionEdge, CompoundGraph> wp = new WeightsFromFile<BioPhysicalEntity, ReactionEdge, CompoundGraph>("/Users/fajourdan/workspace/NetworkComparison/src/main/resources/DistanceMatrix/recon2.v03_ext_noCompartment_noTransport_C-AAM-weights.tab", true);
	        //set weights to edges
	        wp.setWeight(compoundGraph);
	        //remove weights below 0.0
	        int nb = GraphFilter.weightFilter(compoundGraph, 0.0, "<="); System.err.println(nb+" edges removed");
	        //remove edges with NaN weight
	        WeightUtils.removeEdgeWithNaNWeight(compoundGraph);
	        //remove disconnected nodes
	        compoundGraph.removeIsolatedNodes();
	        System.err.println("weights computed.");
	        DefaultWeightPolicy<BioPhysicalEntity, ReactionEdge, CompoundGraph> wpolicy=new DefaultWeightPolicy<BioPhysicalEntity, ReactionEdge, CompoundGraph>();
			wpolicy.setWeight(compoundGraph);
//	        
			HashSet<BioPhysicalEntity> metaboliteToRemove=new HashSet<BioPhysicalEntity>();
			for(BioPhysicalEntity m:metabolites)
			{
				if(!compoundGraph.containsVertex(m))
				{
					metaboliteToRemove.add(m);
					System.err.println("Remove disconnected compound "+ m);
					//System.out.println("Remove disconnected compound "+ m.getName());
				}
			}
			metabolites.removeAll(metaboliteToRemove);
			
//			compoundGraph.removeVertex(compoundGraph.getVertex("M_adp"));
//			compoundGraph.removeVertex(compoundGraph.getVertex("M_atp"));
//			compoundGraph.removeVertex(compoundGraph.getVertex("M_co2"));
//			compoundGraph.removeVertex(compoundGraph.getVertex("M_coa"));
//			compoundGraph.removeVertex(compoundGraph.getVertex("M_ppi"));
//			compoundGraph.removeVertex(compoundGraph.getVertex("M_fadh2"));
//			compoundGraph.removeVertex(compoundGraph.getVertex("M_fad"));
//			compoundGraph.removeVertex(compoundGraph.getVertex("M_h2o"));
//			compoundGraph.removeVertex(compoundGraph.getVertex("M_pi"));
//			compoundGraph.removeVertex(compoundGraph.getVertex("M_nad"));
//			compoundGraph.removeVertex(compoundGraph.getVertex("M_nadh"));
//			compoundGraph.removeVertex(compoundGraph.getVertex("M_nadp"));
//			compoundGraph.removeVertex(compoundGraph.getVertex("M_nadph"));
//			compoundGraph.removeVertex(compoundGraph.getVertex("M_o2"));
//			compoundGraph.removeVertex(compoundGraph.getVertex("M_h"));
//			compoundGraph.removeVertex(compoundGraph.getVertex("M_so4"));
			
			
			HashSet<String> reactionsInPaths=new HashSet<String>();
			
			for(BioPhysicalEntity source : metabolites)
			{
				writer.write(source.getName()+"\t");
				writer.flush();
				int nbMetabolites=metabolites.size();
				int cpt=0;
				
				for(BioPhysicalEntity target:metabolites){
//					if(!compoundGraph.containsVertex(source) || !compoundGraph.containsVertex(target)){
//						if(cpt==nbMetabolites-1){
//							writer.write(Integer.MAX_VALUE);								
//						}
//						else{
//							writer.write(Integer.MAX_VALUE+"\t");
//						}
//					}
//					else{
						if(target.equals(source)){
							if(cpt==nbMetabolites-1){
								writer.write("0");
							}
							else{
								writer.write("0\t");							
							}
						}
						else{
							BioPath<BioPhysicalEntity, ReactionEdge> path=(new ShortestPath<BioPhysicalEntity,ReactionEdge,CompoundGraph>(compoundGraph)).getShortestAsUndirected(source,target);
							
							if(path!=null){
								if(cpt==nbMetabolites-1){
									writer.write(path.getLength()+"");
									//writer.write(path.getLength()*path.getLength()+"");
									//writer.write(path.getWeight()+"");
								}
								else{
									//writer.write(path.getLength()+"\t");
									//writer.write(path.getLength()*path.getLength()+"\t");
									writer.write(path.getWeight()+"\t");
								}
								if (path.getLength()<=3)
								{
									for ( ReactionEdge edge : path.getEdgeList()){
										reactionsInPaths.add(edge.getReaction().getId());
									}
								}
							}
							else{
								if(cpt==nbMetabolites-1){
									writer.write(Integer.MAX_VALUE);								
								}
								else{
									writer.write(Integer.MAX_VALUE+"\t");
								}
							}	
						//}
						
					}
					cpt++;
					writer.flush();					
				}
				writer.write("\n");
				writer.flush();
			}
			

			for(String reaction:reactionsInPaths)
			{
				System.out.println(reaction);
			}
			
			
			
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		

		
	}

}
