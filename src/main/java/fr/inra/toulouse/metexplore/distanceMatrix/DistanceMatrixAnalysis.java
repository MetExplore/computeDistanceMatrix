package fr.inra.toulouse.metexplore.distanceMatrix;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonElement;
import com.google.gson.JsonParser;
import fr.inra.toulouse.metexplore.met4j_core.biodata.*;
import fr.inra.toulouse.metexplore.met4j_core.io.Sbml2BioNetworkLite;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import javax.xml.parsers.ParserConfigurationException;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

/**
 * The Class MeSHAnalysis.
 */
public class DistanceMatrixAnalysis {

    /*
         java -cp "c:/Users/Medday/Documents/Workspace/met4jbinding/target/met4j-binding-1.0-SNAPSHOT-jar-with-dependencies.jar" fr.inra.toulouse.metexplore.Bind2DRank -network "c:/Users/Medday/Documents/results/networkForRank.txt" -edgeWeights "c:/wamp64/www/MetExplore2Conf/atomMapping/testRun-AAM-weights.tab" -metabolites "c:/Users/Medday/Documents/results/seedslist.tab" > c:/Users/Medday/Documents/test/results.json 2> c:/Users/Medday/Documents/test/log & echo $!
    * */
    public String label = "MeSH to Metabolites Analysis";
    public String description = "Computes the Metab2Mesh algorithm to get associated metabolites from a MeSH.";
    public JSONObject networkJson;

    /**
     * metabolic network
     */
    @Option(name = "-fingerprint", usage = "Global network", metaVar = "fingerprintPath", required = true)
    private String fingerprintPath = "fingerprintPath";

    /**
     * metabolic network
     */
    @Option(name = "-network", usage = "Global network", metaVar = "network", required = true)
    private String sbmlPath = "recon2.v03_ext_noCompartment_noTransport.xml";

    /**
     * The objective MeSH
     */
    @Option(name = "-atommapping", usage = "Select MeSH.", metaVar = "recon2.v03_ext_noCompartment_noTransport_C-AAM-weights.tab", required = true)
    private String atomMappingPath = "DistanceMatrix/recon2.v03_ext_noCompartment_noTransport_C-AAM-weights.tab";

    /**
     * The objective MeSH
     */
    @Option(name = "-matrixresult", usage = "Select MeSH.", metaVar = "matrixresult", required = false)
    private String matrixresult;

    /**
     * The objective MeSH
     */
    @Option(name = "-algo", usage = "Select MeSH.", metaVar = "algo", required = false)
    private String algo = "ValidShortest";

    /**
     * The objective MeSH
     */
    @Option(name = "-reactionresult", usage = "Select MeSH.", metaVar = "reactionresult", required = false)
    private String reactionresult;

    /**
     * Constructor
     */
    private DistanceMatrixAnalysis() {
        super();
    }
    public DistanceMatrixAnalysis(String fingerprintPath, String sbmlPath, String atomMappingPath, String algo, String matrixresult, String reactionresult) {
        super();
        this.fingerprintPath = fingerprintPath;
        this.sbmlPath = sbmlPath;
        this.algo = algo;
        this.atomMappingPath = atomMappingPath;
        this.matrixresult = matrixresult;
        this.reactionresult = reactionresult;
    }
    /**
         * The main method.
         *
         * @param args the arguments
         * @throws ParserConfigurationException the parser configuration exception
         * @throws org.xml.sax.SAXException the SAX exception
         * @throws IOException Signals that an I/O exception has occurred.
         */
    public static void main(String[] args) throws ParserConfigurationException, org.xml.sax.SAXException, IOException, IllegalArgumentException, IllegalAccessException {
        DistanceMatrixAnalysis app = new DistanceMatrixAnalysis();

        CmdLineParser parser = new CmdLineParser(app);

        try {
            parser.parseArgument(args);
        } catch (CmdLineException e) {
            e.printStackTrace();
            System.exit(0);
        }
        JSONObject jsonObject = app.run();

        String json = jsonObject.toString();
        System.out.println(jsonObject);
        Gson gson = new GsonBuilder().setPrettyPrinting().create();

        JsonParser jp = new JsonParser();
        JsonElement je = jp.parse(json);
        String prettyJsonString = gson.toJson(je);

        System.out.println(prettyJsonString );
    }

    @SuppressWarnings("unchecked")
    public JSONObject run() throws IllegalArgumentException,
            IllegalAccessException {

        if(this.matrixresult==null)
        {
            System.err.println("Started");
            JSONParser parser = new JSONParser();
            JSONObject networkArg = null;
            try {
                networkArg = (JSONObject) parser.parse(new FileReader(this.sbmlPath));
            } catch (ParseException e) {
                e.printStackTrace();
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            } catch (IOException e) {
                e.printStackTrace();
            }

            this.networkJson = networkArg;
            return this.launchDistanceMatrix();
        }
        return this.launchDistanceMatrixFile();

    }

    protected BioNetwork createBioNetwork() {
        BioNetwork bioNetwork = new BioNetwork();

        JSONArray metabolites = (JSONArray) this.networkJson.get("metabolites");
        for (int i = 0; i < metabolites.size(); ++i) {
            JSONObject metabolite = (JSONObject) metabolites.get(i);
            BioPhysicalEntity bioPhysicalEntity = new BioPhysicalEntity((String) metabolite.get("dbIdentifier"), (String) metabolite.get("name"));
            bioNetwork.addPhysicalEntity(bioPhysicalEntity);
        }

        JSONArray reactions = (JSONArray) this.networkJson.get("reactions");
        for (int i = 0; i < reactions.size(); ++i) {
            JSONObject reaction = (JSONObject) reactions.get(i);

            BioChemicalReaction bioChemicalReaction = new BioChemicalReaction((String)reaction.get("dbIdentifier"), (String)reaction.get("name"));
            bioNetwork.addBiochemicalReaction(bioChemicalReaction);
        }

        HashMap<String, BioChemicalReaction> listReactions =  bioNetwork.getBiochemicalReactionList();
        HashMap<String, BioPhysicalEntity> listMetabolites =  bioNetwork.getPhysicalEntityList();
        JSONArray links = (JSONArray) this.networkJson.get("links");
        for (int i = 0; i < links.size(); ++i) {
            JSONObject link = (JSONObject) links.get(i);
            BioChemicalReaction reaction = listReactions.get(link.get("idReaction"));

            reaction.setReversibility((Boolean) link.get("reversible"));

            if(link.get("interaction").equals("out")){
                BioCompartment compartment = new BioCompartment("comp", "1");
                BioPhysicalEntityParticipant produit = new BioPhysicalEntityParticipant(((Integer) i).toString(), listMetabolites.get(link.get("idMetabolite")), "okok", compartment);

                reaction.addRightParticipant(produit);
            }
            else{
                BioCompartment compartment = new BioCompartment("comp", "1");
                BioPhysicalEntityParticipant substrat = new BioPhysicalEntityParticipant(((Integer) i).toString(), listMetabolites.get(link.get("idMetabolite")), "okok", compartment);

                reaction.addLeftParticipant(substrat);
            }
        }
        return bioNetwork;
    }

    public JSONObject launchDistanceMatrix(){

        DistanceMatrix distanceMatrix = new DistanceMatrix(this.fingerprintPath, this.atomMappingPath, this.algo);
        BioNetwork bioNetwork = this.createBioNetwork();
        HashMap<String, HashMap<String, HashMap>> matrix = distanceMatrix.getMatrixObject(bioNetwork);

        JSONObject json_results_listMetabFromMeSH = new JSONObject(matrix);

        JSONObject jsonObject = new JSONObject();
        jsonObject.put("message", new String("Results are displayed in tables."));
        jsonObject.put("success", "true");
        jsonObject.put("matrix", json_results_listMetabFromMeSH);


        System.err.println("Finished");

        return jsonObject;

    }

    public JSONObject launchDistanceMatrixFile(){

        DistanceMatrix distanceMatrix = new DistanceMatrix(this.fingerprintPath, this.sbmlPath, this.atomMappingPath, this.algo, this.matrixresult, this.reactionresult);
        distanceMatrix.getMatrix();


        JSONObject jsonObject = new JSONObject();
        String rep = new String("Matrix distance is in file: ") + this.matrixresult
                + new String("and path by pairs metabolites are in file: ") + this.reactionresult ;
        jsonObject.put("message", rep);
        jsonObject.put("success", "true");


        System.err.println("Finished");

        return jsonObject;

    }
}
