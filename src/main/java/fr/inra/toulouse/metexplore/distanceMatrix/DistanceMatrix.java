package fr.inra.toulouse.metexplore.distanceMatrix;

import java.io.*;
import java.util.*;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonElement;
import com.google.gson.JsonParser;
import fr.inra.toulouse.metexplore.met4j_core.biodata.BioNetwork;
import fr.inra.toulouse.metexplore.met4j_core.biodata.BioPhysicalEntity;
import fr.inra.toulouse.metexplore.met4j_core.io.Sbml2BioNetworkLite;
import fr.inra.toulouse.metexplore.met4j_graph.computation.algo.KShortestPath;
import fr.inra.toulouse.metexplore.met4j_graph.computation.algo.ShortestPath;
import fr.inra.toulouse.metexplore.met4j_graph.computation.algo.ValidShortestPath;
import fr.inra.toulouse.metexplore.met4j_graph.computation.transform.GraphFilter;
import fr.inra.toulouse.metexplore.met4j_graph.computation.weighting.DefaultWeightPolicy;
import fr.inra.toulouse.metexplore.met4j_graph.computation.weighting.WeightUtils;
import fr.inra.toulouse.metexplore.met4j_graph.computation.weighting.WeightsFromFile;
import fr.inra.toulouse.metexplore.met4j_graph.core.BioPath;
import fr.inra.toulouse.metexplore.met4j_graph.core.compound.CompoundGraph;
import fr.inra.toulouse.metexplore.met4j_graph.core.compound.ReactionEdge;
import fr.inra.toulouse.metexplore.met4j_graph.io.Bionetwork2BioGraph;
import org.json.simple.JSONObject;

public class DistanceMatrix {

	/**
	 * Example args:
	 * E:\Recherche\RESEARCH-DB-Networks\Homo_sapiens\recon2.v03_ext_noCompartment_noTransport.xml E:\workspace\NetworkComparison\src\main\resources\DistanceMatrix\Sonia\fingerprint.txt E:\workspace\NetworkComparison\src\main\resources\DistanceMatrix\Sonia\distanceMatrix.txt
	 * @param args
	 */

	private String fingerprint;
	private String sbmlPath;
	private String atomMappingFileName;
	private String matrixPath;
	private String reactionsPath;
	private String algo;

	private BioNetwork bionetwork;

	public DistanceMatrix(String fingerprintPath, String atomMappingFileName, String algo) {

		this.fingerprint = fingerprintPath;
		this.atomMappingFileName = atomMappingFileName;
		this.algo = algo;
		this.matrixPath = null;
		this.reactionsPath = null;
	}

	public DistanceMatrix(String fingerprintPath, String sbmlPath, String atomMappingFileName, String algo, String matrixPath, String reactionsPath) {
		this.fingerprint = fingerprintPath;
		this.sbmlPath = sbmlPath;
		this.algo = algo;
		this.atomMappingFileName = atomMappingFileName;
		this.matrixPath = matrixPath;
		this.reactionsPath = reactionsPath;
	}

	public void getMatrix(){

		System.err.println("----GO");

		/**
		 * Graph creation
		 */
		//import the network
		Sbml2BioNetworkLite importer=new Sbml2BioNetworkLite(this.sbmlPath,true);
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
			BufferedReader reader=new BufferedReader(new FileReader(this.fingerprint));
			try {
				String line=reader.readLine();
				while(line!=null)
				{
					String[] elements=line.split("\t");
					BioPhysicalEntity metabolite=net.getBioPhysicalEntityById(elements[0]);
					metabolites.add(metabolite);
					if (metabolite==null)
						System.err.println(elements[0]+" not found ");
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
			BufferedWriter writer=new BufferedWriter(new FileWriter(this.matrixPath));
			BufferedWriter writerReactions=new BufferedWriter(new FileWriter(this.reactionsPath));

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
			WeightsFromFile<BioPhysicalEntity, ReactionEdge, CompoundGraph> wp = new WeightsFromFile<BioPhysicalEntity, ReactionEdge, CompoundGraph>(this.atomMappingFileName, true);
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

			HashSet<String> reactionsInPaths=new HashSet<String>();
			HashMap<String, HashSet> paths = new HashMap<String, HashSet>();

            System.err.println("Start computing shortest paths.");

			HashMap<String, Integer> dist = new HashMap<String, Integer>();
			for(BioPhysicalEntity source : metabolites)
			{
				int nbMetabolites=metabolites.size();
				int cpt=0;
				
				for(BioPhysicalEntity target:metabolites){

					int compare = source.getId().compareTo(target.getId());
					String pathName;
					if (compare < 0) {
						pathName = source.getId() + target.getId();
					} else {
						pathName = target.getId() + source.getId();
					}
					if (!dist.containsKey(pathName)){

						if(!compoundGraph.containsVertex(source) || !compoundGraph.containsVertex(target)){
								dist.put(pathName, Integer.MAX_VALUE);
						}
						else{
							if(target.equals(source)){
								dist.put(pathName, 0);
							}
							else {
								BioPath<BioPhysicalEntity, ReactionEdge> pathsimple = null;
								BioPath<BioPhysicalEntity, ReactionEdge> pathreverse = null;

								if (this.algo.equals("ShortestAsUndirected"))
								{
									pathsimple = (new ShortestPath<BioPhysicalEntity, ReactionEdge, CompoundGraph>(compoundGraph)).getShortestAsUndirected(source, target);
									pathsimple = (new ShortestPath<BioPhysicalEntity, ReactionEdge, CompoundGraph>(compoundGraph)).getShortestAsUndirected(target, source);

								}
								else{
									pathsimple = (new ValidShortestPath(compoundGraph)).getValidShortest(source, target, 15);
									pathreverse = (new ValidShortestPath(compoundGraph)).getValidShortest(target, source, 15);
								}
								BioPath<BioPhysicalEntity, ReactionEdge> path = null;

								if (pathsimple != null || pathreverse != null){
									if(pathsimple == null){
										path = pathreverse;
									}
									else {
										if(pathreverse == null){
											path = pathsimple;
										}
										else
										{
											if (pathreverse.getLength()<pathsimple.getLength()){
												path = pathreverse;
											}
											else{
												path = pathsimple;
											}
										}
									}
									dist.put(pathName, path.getLength());


									if (!paths.keySet().contains(pathName)) {
										HashSet<String> reactionsInPath = new HashSet<String>();
										for (ReactionEdge edge : path.getEdgeList()) {
											reactionsInPaths.add(edge.getReaction().getId());
											reactionsInPath.add(edge.getReaction().getId());
										}
										paths.put(pathName, reactionsInPath);
									}
								}
								else {
									System.err.println("No valid path found between "+source.getName()+" and "+target.getName()+".");
									dist.put(pathName, Integer.MAX_VALUE);
								}
							}















//							BioPath<BioPhysicalEntity, ReactionEdge> path=(new ShortestPath<BioPhysicalEntity,ReactionEdge,CompoundGraph>(compoundGraph)).getShortestAsUndirected(source,target);
//							List<BioPath<BioPhysicalEntity, ReactionEdge>> kpaths = (new KShortestPath<BioPhysicalEntity,ReactionEdge,CompoundGraph>(compoundGraph)).getKShortest(source,target, 10);
//
//							boolean found = false;
//                            if(kpaths.size()>0){
//                                Iterator<BioPath<BioPhysicalEntity, ReactionEdge>> iterator = kpaths.iterator();
//                                while (iterator.hasNext() && !found){
//                                    BioPath<BioPhysicalEntity, ReactionEdge> path = iterator.next();
//                                    System.out.println(ValidShortestPath.isValid(path));
//                                    System.out.println(path);
//                                    if(ValidShortestPath.isValid(path)){
//                                        found = true;
//                                        if(cpt==nbMetabolites-1){
//                                            writer.write(path.getLength()+"");
//                                            //writer.write(path.getLength()*path.getLength()+"");
//                                            //writer.write(path.getWeight()+"");
//                                        }
//                                        else{
//                                            //writer.write(path.getLength()+"\t");
//                                            //writer.write(path.getLength()*path.getLength()+"\t");
//                                            writer.write(path.getLength()+"\t");
//                                        }
//
//                                        int compare = source.getId().compareTo(target.getId());
//                                        String pathName;
//                                        if (compare < 0){
//                                            pathName = source.getId()+target.getId();
//                                        }
//                                        else
//                                        {
//                                            pathName = target.getId()+source.getId();
//                                        }
//                                        if(!paths.keySet().contains(pathName)){
//                                            HashSet<String> reactionsInPath = new HashSet<String>();
//                                            for ( ReactionEdge edge : path.getEdgeList()){
//                                                reactionsInPaths.add(edge.getReaction().getId());
//                                                reactionsInPath.add(edge.getReaction().getId());
//                                            }
//                                            paths.put(pathName, reactionsInPath);
//                                        }
//                                    }
//                                }
//                                if(!found){
//                                    if(cpt==nbMetabolites-1){
//                                        writer.write(Integer.MAX_VALUE);
//                                    }
//                                    else{
//                                        writer.write(Integer.MAX_VALUE+"\t");
//                                    }
//                                    System.out.println(kpaths);
//                                    System.err.println("No valid path found between "+source.getName()+" and "+target.getName()+". Choose highest k for k-shortest path.");
//                                }
//                            }
//                            else{
//                                if(cpt==nbMetabolites-1){
//                                    writer.write(Integer.MAX_VALUE);
//                                }
//                                else{
//                                    writer.write(Integer.MAX_VALUE+"\t");
//                                }
//                            }

							//}

						}

					}

				}
			}
			for(BioPhysicalEntity source : metabolites)
			{

				writer.write(source.getName()+"\t");
				writer.flush();
				int nbMetabolites=metabolites.size();
				int cpt=0;

				for(BioPhysicalEntity target: metabolites) {

					int compare = source.getId().compareTo(target.getId());
					String pathName;
					if (compare < 0) {
						pathName = source.getId() + target.getId();
					} else {
						pathName = target.getId() + source.getId();
					}

					if(cpt==nbMetabolites-1){
						writer.write(dist.get(pathName)+"\n");
					}
					else{
						writer.write(dist.get(pathName)+"\t");
					}

					cpt++;
				}
				writer.flush();
			}

            Gson gson = new GsonBuilder().setPrettyPrinting().create();

            JsonParser jp = new JsonParser();
            JsonElement je = jp.parse(new JSONObject(paths).toString());
            String prettyJsonString = gson.toJson(je);

            writerReactions.write(prettyJsonString);
            writerReactions.flush();
			
//
//			for(String reaction:reactionsInPaths)
//			{
//				System.out.println(reaction);
//			}
//
			
			
			writer.close();
			writerReactions.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public HashMap<String, HashMap<String, HashMap>> getMatrixObject(BioNetwork net){

		System.err.println("----GO");
        /**
         * Graph creation
         */
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
            BufferedReader reader=new BufferedReader(new FileReader(this.fingerprint));
            try {
                String line=reader.readLine();
                while(line!=null)
                {
                    String[] elements=line.split("\t");
                    BioPhysicalEntity metabolite=net.getBioPhysicalEntityById(elements[0]);
                    metabolites.add(metabolite);
                    if (metabolite==null)
                        System.err.println(elements[0]+" not found ");
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
		 * Main part of the function that will create a distance matrix
		 */
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
			WeightsFromFile<BioPhysicalEntity, ReactionEdge, CompoundGraph> wp = new WeightsFromFile<BioPhysicalEntity, ReactionEdge, CompoundGraph>(this.atomMappingFileName, true);
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

			HashMap<String, HashMap<String, HashMap>> results = new HashMap<String, HashMap<String, HashMap>>();

			HashSet<String>  reactionsInPaths = new HashSet<String> ();
			for(BioPhysicalEntity source : metabolites)
			{
				if(!results.containsKey(source.getName()))
					results.put(source.getName(), new HashMap());

				int nbMetabolites=metabolites.size();
				int cpt=0;

				for(BioPhysicalEntity target:metabolites){
                    if(!results.get(source.getName()).containsKey(target.getName()))
                        results.get(source.getName()).put(target.getName(), new HashMap());
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
                        results.get(source.getName()).get(target.getName()).put("score", 0);
					}
					else{
						BioPath<BioPhysicalEntity, ReactionEdge> path=(new ShortestPath<BioPhysicalEntity,ReactionEdge,CompoundGraph>(compoundGraph)).getShortestAsUndirected(source,target);

						if(path!=null){
                            results.get(source.getName()).get(target.getName()).put("score", path.getLength());

							HashSet<String> reactionsArray = new HashSet<String>();
							for ( ReactionEdge edge : path.getEdgeList()){
                                reactionsArray.add(edge.getReaction().getId());
								reactionsInPaths.add(edge.getReaction().getId());
							}
                            results.get(source.getName()).get(target.getName()).put("reactions", reactionsArray);
						}
						else{
                            results.get(source.getName()).get(target.getName()).put("score", Integer.MAX_VALUE);
						}
						//}

					}
					cpt++;
				}
			}


			for(String reaction:reactionsInPaths)
			{
				System.err.println(reaction);
			}

		return results;
	}

	private String getFile(String fileName) {

		StringBuilder result = new StringBuilder("");

		//Get file from resources folder
		ClassLoader classLoader = getClass().getClassLoader();
		File file = new File(classLoader.getResource(fileName).getFile());
		return file.toString();
	}
}
