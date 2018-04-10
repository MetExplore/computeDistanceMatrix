package fr.inra.toulouse.metexplore;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import org.kohsuke.args4j.Option;

import fr.inra.toulouse.metexplore.met4j_core.biodata.BioChemicalReaction;
import fr.inra.toulouse.metexplore.met4j_core.biodata.BioNetwork;
import fr.inra.toulouse.metexplore.met4j_core.biodata.BioPhysicalEntity;
import fr.inra.toulouse.metexplore.met4j_core.biodata.BioPhysicalEntityParticipant;
import fr.inra.toulouse.metexplore.met4j_core.io.Sbml2BioNetworkLite;
import fr.inra.toulouse.metexplore.met4j_graph.computation.algo.SteinerTreeApprox;
import fr.inra.toulouse.metexplore.met4j_graph.computation.analysis.ConnectedComponents;
import fr.inra.toulouse.metexplore.met4j_graph.computation.analysis.GraphMeasure;
import fr.inra.toulouse.metexplore.met4j_graph.core.compound.CompoundGraph;
import fr.inra.toulouse.metexplore.met4j_graph.core.reaction.CompoundEdge;
import fr.inra.toulouse.metexplore.met4j_graph.core.reaction.ReactionGraph;
import fr.inra.toulouse.metexplore.met4j_graph.io.Bionetwork2BioGraph;
import fr.inra.toulouse.metexplore.met4j_graph.io.ExportGraph;
import fr.inra.toulouse.metexplore.met4j_toolbox.AbstractApplication;



/**
 * @author Fabien Jourdan
 * Function to steiner tree connecting reactions moving from one status to another
 * Given a metabolic network and information on status change of reactions between two conditions
 * Reactions which move from an initial status (e.g. inactivated) to a final status (e.g. activated) 
 *  will be grouped in a tree (Steiner tree)
 */
public class ExtractSteinerTreeSwitchingReactions extends AbstractApplication{
	@Override
	public String getDescription() {
		// TODO Auto-generated method stub
		String description="Function to steiner tree connecting reactions moving from one status to another\n"
				+ "Given a metabolic network and information on status change of reactions between two conditions \n"
				+ "Reactions which move from an initial status (e.g. inactivated) to a final status (e.g. activated)\n"
				+ "will be grouped in a tree (Steiner tree)";
		return description;
	}

	@Option(name="-s", usage="[sbmlFile] Sbml file -- Required")
	private String sbmlFile;
	public String getSbmlFile() {
		return sbmlFile;
	}
	
	@Option(name="-st", usage="[statusFile] File containing status with per columns: ReactionId -1 (inactivated) 0 (both) 1 (activated) -- Required")
	private String statusFile;
	public String getStatusFile() {
		return statusFile;
	}

	@Option(name="-del", usage="[removeInactiveReactions] If true, remove reactions not active (-1) at both stages. -- Optinoal")
	private Boolean removeInactiveReactions=true;
	public Boolean getRemoveInactiveReactions() {
		return removeInactiveReactions;
	}

	
	@Option(name="-stsep", usage="[statusSeparator] Separator used in status file -- Optional")
	private String statusSeparator=";";
	public String getStatusSeparator() {
		return statusSeparator;
	}	
	
	
	@Option(name="-oC", usage="[connectedComponentsOutputFile] Output file for connected components-- Optional")
	private String connectedComponentsOutputFile="connectedComponents.txt";
	public String getConnectedComponentsOutputFile() {
		return connectedComponentsOutputFile;
	}
	
	@Option(name="-oSt", usage="[steinerOutputFile] Output file for Steiner Tree-- Optional")
	private String steinerOutputFile="steinerTree.txt";
	public String getSteinerOutputFile() {
		return steinerOutputFile;
	}
	
	@Option(name="-from", usage="[initialStatus] Initial reaction status -- Optional")
	private String initialStatus="-1";
	public String getInitialStatus() {
		return initialStatus;
	}	
	
	@Option(name="-to", usage="[initialStatus] Final reaction stats -- Optional")
	private String finalStatus="1";
	public String getFinalStatus() {
		return finalStatus;
	}	
	
	@Option(name="-sc", usage="[optional] File containing side compound identifiers in a single column -- Optional")
	private String sideCompoundsFile;
	public String getSideCompoundsFile() {
		return sideCompoundsFile;
	}


	public static void main(String[] args) {
	//create the object to get the arguments handled
		ExtractSteinerTreeSwitchingReactions extractor=new ExtractSteinerTreeSwitchingReactions();
		extractor.parseArguments(args);
	//import the network
		Sbml2BioNetworkLite importer=new Sbml2BioNetworkLite(extractor.getSbmlFile(),false);
		System.out.println("-- SBML imported");
	//Create the bionetwork
		BioNetwork network=importer.getBioNetwork();
		System.out.println("-- BioNetwork created");
	//Remove side compounds from the bionetwork if a file is provide
		if(extractor.getSideCompoundsFile()!=null){
			removeSideCompounds(network, extractor.getSideCompoundsFile());	
			System.out.println("-- Side compounds removed");
		}
	//Get reactions activities, create a map with for each node an array for its different status
		HashMap<BioChemicalReaction,String[]> reactionStatusMap=new HashMap<BioChemicalReaction, String[]>();
		try {
			BufferedReader reader=new BufferedReader(new FileReader(extractor.getStatusFile()));
			try {
				String line=reader.readLine();
				while(line!=null)
				{
					String[] elements=line.split(extractor.getStatusSeparator());
					BioChemicalReaction reaction=network.getBiochemicalReactionList().get(elements[0]);
					reactionStatusMap.put(reaction,Arrays.copyOfRange(elements, 1, 3));
					line=reader.readLine();
				}
			} catch (IOException e) {
				e.printStackTrace();
			}
		} catch (FileNotFoundException e) {
			System.err.println("Please provide status file");
			e.printStackTrace();
		}
		
	//Remove reactions :
	//		- that are not in the input list of reactions with status
	//		- that are not in active in both conditions (-1) if it is asked by user
	//		- reactions that are not active in the second condition
		ArrayList<BioChemicalReaction> reactionsToRemove=new ArrayList<BioChemicalReaction>();
		for(BioChemicalReaction reaction:network.getBiochemicalReactionList().values()){
			//if there is no activity status information on a reaction: remove it
			if (reactionStatusMap.get(reaction)==null){
				//System.out.println(reactionStatusMap.get(reaction));
				reactionsToRemove.add(reaction);
			}
			else{
				if(extractor.getRemoveInactiveReactions())
				{
					if(
							((reactionStatusMap.get(reaction)[0].equals("-1"))&&
							(reactionStatusMap.get(reaction)[1].equals("-1")))
							){
//						System.err.println("remove ..."+reaction.getId());
						reactionsToRemove.add(reaction);					
					}
					if(reactionStatusMap.get(reaction)[1].equals("-1")){
						reactionsToRemove.add(reaction);
					}
					
				}
			}
		}
		
		System.err.println("nb of reactions "+network.getBiochemicalReactionList().values().size());
		
		for(BioChemicalReaction reaction : reactionsToRemove)
		{
			network.removeBioChemicalReaction(reaction.getId());
		}		

		System.err.println("nb of reactions "+network.getBiochemicalReactionList().values().size());
		
		
		
	//Remove reactions that are not moving from initial status to final status
//		ArrayList<BioChemicalReaction> reactionsToRemove=new ArrayList<BioChemicalReaction>();
//		for(BioChemicalReaction reaction:network.getBiochemicalReactionList().values()){
//			//if there is no activity status information on a reaction: remove it
//			if (reactionStatusMap.get(reaction)==null){
//				//System.out.println(reactionStatusMap.get(reaction));
//				reactionsToRemove.add(reaction);
//			}
//			else{
//				if(
//						!((reactionStatusMap.get(reaction)[0].equals(extractor.getInitialStatus()))&&
//						(reactionStatusMap.get(reaction)[1].equals(extractor.getFinalStatus())))
//						){
//					reactionsToRemove.add(reaction);					
//				}
//			}
//		}
//		for(BioChemicalReaction reaction : reactionsToRemove)
//		{
//			network.removeBioChemicalReaction(reaction.getId());
//		}
	//Build reaction graph from remaining network
		ReactionGraph reactionGraph=buildReactionGraph(network);
		System.out.println("-- Reaction graph created");
		System.err.println("Number of nodes in reaction graph "+reactionGraph.vertexSet().size());
		System.err.println("Number of edges in reaction graph "+reactionGraph.edgeSet().size());
		
	//Get list of reactions with the right status change and the other ones
		HashSet<BioChemicalReaction> targets=new HashSet<BioChemicalReaction>();
		HashSet<BioChemicalReaction> otherStatus=new HashSet<BioChemicalReaction>();
		for(BioChemicalReaction reaction : reactionGraph.vertexSet())
		{
			if(
					((reactionStatusMap.get(reaction)[0].equals(extractor.getInitialStatus()))&&
					(reactionStatusMap.get(reaction)[1].equals(extractor.getFinalStatus())))
					){
				targets.add(reaction);					
			}
			else{
				otherStatus.add(reaction);
			}
		}
		System.out.println("-- Switching reactions identified");
		//Create Connected components file
		System.out.println("Number of reactions to remove "+otherStatus.size());
		for(BioChemicalReaction reaction:otherStatus)
		{
			network.removeBioChemicalReaction(reaction.getId());
		}
		ReactionGraph reactionGraphForConnectedC=buildReactionGraph(network);
		System.out.println("-- Filtered reaction graph created");
		//reactionGraphForConnectedC.removeAllVertices(otherStatus);
		System.err.println("Number of nodes in filtered reaction graph "+reactionGraphForConnectedC.vertexSet().size());
		System.err.println("Number of edges in filtered reaction graph "+reactionGraphForConnectedC.edgeSet().size());
		computeAndPrintConnectedComponentResults(reactionGraphForConnectedC,extractor.getConnectedComponentsOutputFile());
		
		
		System.err.println("Number of seeds "+targets.size());
		System.err.println("Number of nodes in reaction graph "+reactionGraph.vertexSet().size());
		System.err.println("Number of edges in reaction graph "+reactionGraph.edgeSet().size());
	//compute Steiner tree
		SteinerTreeApprox<BioChemicalReaction, CompoundEdge, ReactionGraph> steinerCompute=new SteinerTreeApprox<BioChemicalReaction, CompoundEdge, ReactionGraph>(reactionGraph);
		ArrayList<CompoundEdge> steinerEdgeSet=new ArrayList<CompoundEdge>(steinerCompute.getSteinerTreeList(targets,false));
		System.out.println(steinerCompute.getSteinerTreeList(targets,false));
		System.out.println("Number of connected components "+GraphMeasure.getConnectedCompenent(reactionGraph).size());
		printSteinerResult(reactionGraph, extractor.getSteinerOutputFile(), steinerEdgeSet,targets);
		
		
	}
	
	static void printSteinerResult(ReactionGraph graph,String outputFilePath,ArrayList<CompoundEdge> compoundEdges,HashSet<BioChemicalReaction> targets)
	{
		HashSet<BioChemicalReaction> reactionsInTree=new HashSet<BioChemicalReaction>();
		for(CompoundEdge edge : compoundEdges)
		{
			reactionsInTree.add(edge.getV1());
			reactionsInTree.add(edge.getV2());
		}
		System.err.println("Nob of reactions in tree "+reactionsInTree.size());
		try {
			BufferedWriter writer=new BufferedWriter(new FileWriter(outputFilePath));
			for(BioChemicalReaction reaction : graph.vertexSet())
			{
				if(targets.contains(reaction)){
					writer.write(reaction.getId()+"\t"+2+"\n");					
				}
				else{
					if(reactionsInTree.contains(reaction)){
						writer.write(reaction.getId()+"\t"+1+"\n");	
					}
					else
					{
						writer.write(reaction.getId()+"\t"+0+"\n");						
					}
				}
					writer.flush();
			}
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	
	
	static ReactionGraph buildReactionGraph(BioNetwork network)
	{
		ReactionGraph reactionGraph = new ReactionGraph();
		HashSet<BioChemicalReaction> exchange = new HashSet<BioChemicalReaction>();
		for(BioChemicalReaction r : network.getBiochemicalReactionList().values()){
			if(!r.isExchangeReaction()){
				reactionGraph.addVertex(r);
			}else{
				exchange.add(r);
			}
		}
		//for each reaction r1 in N
		int cpt=0;
		for(BioChemicalReaction r1 : network.getBiochemicalReactionList().values())
		{
			System.out.println(" "+cpt++);
			//for each reaction r2 in N
			for(BioChemicalReaction r2 : network.getBiochemicalReactionList().values())
			{
				//if the two reactions are different and not exchange ones
				if((r1!=r2)&&(!exchange.contains(r1))&&(!exchange.contains(r2))){
					//for each mR1 in right participant of r1
					for(BioPhysicalEntityParticipant mR1 : r1.rightParticipantList.values()){
						//for each mL2 in left participant of r2
						for(BioPhysicalEntityParticipant mL2 : r2.leftParticipantList.values()){
							if(mR1.getPhysicalEntity().getId().equals(mL2.getPhysicalEntity().getId())){
								reactionGraph.addVertex(r1);
								reactionGraph.addVertex(r2);
								reactionGraph.addEdge(r1, r2, new CompoundEdge(r1,r2,mR1.getPhysicalEntity()));
							}
						}
						//if r2 is reversible
						if(r2.isReversible())
						{
							//for each mR2 in right participant of r2
							for(BioPhysicalEntityParticipant mR2 : r2.rightParticipantList.values()){
								if(mR1.getPhysicalEntity().getId().equals(mR2.getPhysicalEntity().getId())){
									reactionGraph.addVertex(r1);
									reactionGraph.addVertex(r2);
									reactionGraph.addEdge(r1, r2, new CompoundEdge(r1,r2,mR1.getPhysicalEntity()));
								}
							}
						}						
					}
					if(r1.isReversible())
					{
						for(BioPhysicalEntityParticipant mL1 : r1.leftParticipantList.values()){
							//for each mL2 in left participant of r2
							for(BioPhysicalEntityParticipant mL2 : r2.leftParticipantList.values()){
								if(mL1.getPhysicalEntity().getId().equals(mL2.getPhysicalEntity().getId())){
									reactionGraph.addVertex(r1);
									reactionGraph.addVertex(r2);
									reactionGraph.addEdge(r1, r2, new CompoundEdge(r1,r2,mL1.getPhysicalEntity()));
								}
							}
							//if r2 is reversible
							if(r2.isReversible())
							{
								//for each mR2 in right participant of r2
								for(BioPhysicalEntityParticipant mR2 : r2.rightParticipantList.values()){
									if(mL1.getPhysicalEntity().getId().equals(mR2.getPhysicalEntity().getId())){
										reactionGraph.addVertex(r1);
										reactionGraph.addVertex(r2);
										reactionGraph.addEdge(r1, r2, new CompoundEdge(r1,r2,mL1.getPhysicalEntity()));
									}
								}
							}						
						}						
					}
				}
			}
		}
		return reactionGraph;
	}
	
	static void removeSideCompounds(BioNetwork network,String sideCompoundsFile)
	{
		try {
			BufferedReader reader=new BufferedReader(new FileReader(sideCompoundsFile));
			try {
				String line=reader.readLine();
				while(line!=null)
				{
					String[] elements=line.split("\t");
					network.removeCompound(elements[0]);
					line=reader.readLine();
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	static void computeAndPrintConnectedComponentResults(ReactionGraph graph,String outputFilePath)
	{
		try {
			BufferedWriter writer=new BufferedWriter(new FileWriter(outputFilePath));
			int componentNB=0;
			for(Set<BioChemicalReaction> component : GraphMeasure.getConnectedCompenent(graph))
			{
				for(BioChemicalReaction reaction: component)
				{
					writer.write(reaction.getId()+"\t"+componentNB+"\n");
					writer.flush();
				}
				componentNB++;			
			}
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	static void test(){
		// TODO Auto-generated method stub
		ReactionGraph testGraph=new ReactionGraph();
		BioChemicalReaction A=new BioChemicalReaction("A");
		BioChemicalReaction B=new BioChemicalReaction("B");
		BioChemicalReaction C=new BioChemicalReaction("C");
		BioChemicalReaction D=new BioChemicalReaction("D");
		BioChemicalReaction E=new BioChemicalReaction("E");
		BioChemicalReaction F=new BioChemicalReaction("F");
		BioChemicalReaction G=new BioChemicalReaction("G");
		BioChemicalReaction H=new BioChemicalReaction("H");
		BioChemicalReaction I=new BioChemicalReaction("I");
		testGraph.addVertex(A);
		testGraph.addVertex(B);
		testGraph.addVertex(C);
		testGraph.addVertex(D);
		testGraph.addVertex(E);
		testGraph.addVertex(F);
		testGraph.addVertex(G);
		testGraph.addVertex(H);
		testGraph.addVertex(I);
		BioPhysicalEntity c1=new BioPhysicalEntity("c1");
		BioPhysicalEntity c2=new BioPhysicalEntity("c2");
		BioPhysicalEntity c3=new BioPhysicalEntity("c3");
		BioPhysicalEntity c4=new BioPhysicalEntity("c4");
		BioPhysicalEntity c5=new BioPhysicalEntity("c5");
		BioPhysicalEntity c6=new BioPhysicalEntity("c6");
		BioPhysicalEntity c7=new BioPhysicalEntity("c7");
		BioPhysicalEntity c8=new BioPhysicalEntity("c8");
		BioPhysicalEntity c9=new BioPhysicalEntity("c9");
		CompoundEdge AB=new CompoundEdge(A,B, c2);
		CompoundEdge BC=new CompoundEdge(B,C, c3);
		CompoundEdge BD=new CompoundEdge(B,D, c3);
		CompoundEdge CE=new CompoundEdge(C,E, c4);
		CompoundEdge CF=new CompoundEdge(C,F, c4);
		CompoundEdge DG=new CompoundEdge(D,G, c5);
		CompoundEdge DH=new CompoundEdge(D,H, c5);
		CompoundEdge FI=new CompoundEdge(F,I, c7);
		CompoundEdge GI=new CompoundEdge(G,I,c7);
		testGraph.addEdge(A,B,AB);
		testGraph.addEdge(B,C,BC);
		testGraph.addEdge(B,D,BD);
		testGraph.addEdge(C,E,CE);
		testGraph.addEdge(C,F,CF);
		testGraph.addEdge(D,G,DG);
		testGraph.addEdge(D,H,DH);
		testGraph.addEdge(F,I,FI);
		testGraph.addEdge(G,I,GI);
		System.out.println("Graph number of edges : "+testGraph.edgeSet().size());
		System.out.println("Graph number of nodes : "+testGraph.vertexSet().size());
		
		HashMap<BioChemicalReaction,String[]> reactionStatusMap=new HashMap<BioChemicalReaction, String[]>();
		//create a map with for each node an array for its different status
		importStatus("E:\\workspace\\NetworkComparison\\src\\test\\resources\\status.txt","\t",testGraph,reactionStatusMap);
		
		//Remove all nodes which doesn't go from 1 to -1
		ArrayList<BioChemicalReaction> reactionsToRemove=new ArrayList<BioChemicalReaction>();
		for (BioChemicalReaction reactionNode : testGraph.vertexSet())
		{
			if(
					!((reactionStatusMap.get(reactionNode)[0].equals("1")) &&
					(reactionStatusMap.get(reactionNode)[1].equals("-1")))
					){
				reactionsToRemove.add(reactionNode);
			}
		}
		for(BioChemicalReaction r : reactionsToRemove)
		{
			testGraph.removeVertex(r);
		}
		
		System.out.println("Graph number of edges : "+testGraph.edgeSet().size());
		System.out.println("Graph number of nodes : "+testGraph.vertexSet().size());
		System.out.println("Number of connected components "+GraphMeasure.getConnectedCompenent(testGraph).size());
		
		computeAndPrintConnectedComponentResults(testGraph,"E:\\workspace\\NetworkComparison\\src\\test\\resources\\result.txt");
	}
	
	static void importStatus(String filePath, String separator,ReactionGraph reactionGraph,HashMap<BioChemicalReaction,String[]> reactionStatusMap)
	{
		try {
			BufferedReader reader=new BufferedReader(new FileReader(filePath));
			try {
				String line=reader.readLine();
				while(line!=null)
				{
					String[] elements=line.split("\t");
					BioChemicalReaction reaction=reactionGraph.getVertex(elements[0]);
					reactionStatusMap.put(reaction,Arrays.copyOfRange(elements, 1, 3));
					line=reader.readLine();
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}



}
