package fr.inra.toulouse.metexplore.distanceMatrix;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonElement;
import com.google.gson.JsonParser;
import com.sun.scenario.effect.impl.sw.sse.SSEBlend_SRC_OUTPeer;
import fr.inra.toulouse.metexplore.met4j_core.biodata.*;
import fr.inra.toulouse.metexplore.met4j_core.io.Sbml2BioNetworkLite;
import org.apache.maven.model.Model;
import org.apache.maven.model.io.xpp3.MavenXpp3Reader;
import org.codehaus.plexus.util.xml.pull.XmlPullParserException;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import javax.xml.parsers.ParserConfigurationException;
import java.io.*;
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
    @Option(name = "-h", usage = "Help", metaVar = "Help", required = false)
    private boolean help = false;

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
     * Arg to use the class in MetExplore
     *
     */
    @Option(name = "-useMetExploreOutput", usage = "useMetExploreOutput", metaVar = "useMetExploreOutput", required = false)
    private boolean useMetExploreOutput = false;



    /**
     * Constructor
     */
    private DistanceMatrixAnalysis() {
        super();
    }
    public DistanceMatrixAnalysis(String fingerprintPath, String sbmlPath, String atomMappingPath, String algo, String matrixresult, String reactionresult, String useMetExploreOutput, String help) {
        super();
        this.fingerprintPath = fingerprintPath;
        this.sbmlPath = sbmlPath;
        this.algo = algo;
        this.atomMappingPath = atomMappingPath;
        this.matrixresult = matrixresult;
        this.reactionresult = reactionresult;

        this.useMetExploreOutput = (useMetExploreOutput.toLowerCase().equals("true")) ? true : false ;
        this.help = (help.toLowerCase().equals("true")) ? true : false ;
    }
    /**
         * The main method.
         *
         * @param args the arguments
         * @throws ParserConfigurationException the parser configuration exception
         * @throws org.xml.sax.SAXException the SAX exception
         * @throws IOException Signals that an I/O exception has occurred.
         */
    public static void main(String[] args) throws ParserConfigurationException, org.xml.sax.SAXException, IOException, IllegalArgumentException, IllegalAccessException, XmlPullParserException {
        DistanceMatrixAnalysis app = new DistanceMatrixAnalysis();

        CmdLineParser parser = new CmdLineParser(app);

        try {
            parser.parseArgument(args);
        } catch (CmdLineException e) {
            e.printStackTrace();
            System.exit(0);
        }
        app.run();


    }

    @SuppressWarnings("unchecked")
    public void run() throws IllegalArgumentException,
            IllegalAccessException, IOException, XmlPullParserException {

        if(this.help){
            MavenXpp3Reader reader = new MavenXpp3Reader();
            Model model = reader.read(new FileReader("pom.xml"));
            System.out.println( "                           `-:+syhdmNNNNNNNNmdhys+:-`                           ");
            System.out.println( "                      `-+ydNNmhso+/:-......--:/osydNNdy+-`                      ");
            System.out.println( "                   -ohNNds/-`                       ./ohmNho-                   ");
            System.out.println( "                -omNds:`          `.-::::::--`          `-ohNms-                ");
            System.out.println( "             `+dNdo-       `-/oydmNNNNNNNNNNy   /s`         `+hNd+`             ");
            System.out.println( "           .omNy:`      -+hmNNNNNNNNNNNNNNNNm/-:dNy.   -`      .ommo.           ");
            System.out.println( "         `omNs.      :smNNNNNNNNNNNNNNNNNNNNNhsshNNNdhdNNh/`     `+mmo`         ");
            System.out.println( "        /mNy.     `+dNNNNNNNNNNNNNNNNNNNNNNN:`  `:NNNs`.yNNms-     `omm/        ");
            System.out.println( "      `yNd:     `+mNNNNNNNNNNNNNNNNNNNNNNNNm      mNNh/-sdNNNNy-     .yNy`      ");
            System.out.println( "     -dNy`     /dNNNNNNNNNNNNNNNNNNNNNNNNNNNh/--/hdy+:````:mNNNNs`     +Nd-     ");
            System.out.println( "    :mN+     `yNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNms/.`        -NNNNNm:     -mm:    ");
            System.out.println( "   :mN+     -dNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN/             /NNNNNNo     -mN:   ");
            System.out.println( "  .mN+     -mNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNm               oNNNNNNs     -mm.  ");
            System.out.println( " `hNy     .mNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNm       /        yNNNNNNo     /Nd` ");
            System.out.println( " +Nm`    `dNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN`      d:       `hNNNNNN:     yN+ ");
            System.out.println( "`mN+     +NNNNNNNNNNNNmddmNNNNNNNNNNNNNNNNNNNN.      hm.       .dNNNNNd     .Nm`");
            System.out.println( "/NN`     mNNNNNNNmhs/-.``.:yNNNNNNNNNNNNNNNNNN-      sNd`       -mNNNNN/     yN/");
            System.out.println( "yNh     :NNNNNy+-.          :hNNNNNNNNNNNNNNNN:      +NNh`       :NNNNNy     /Ny");
            System.out.println( "dNo     +NNNN/               `/dNNNNNNNNNNNNNN/      /NNNs        +NNNNm     -Nm");
            System.out.println( "NN+     sNNNN+        ``       `+mNNNNNNNNNNNN+      -NNNN+        sNNNN     `NN");
            System.out.println( "NN+     oNNNNN-       `o:        .omNNNNNNNNNNo      .NNNNN/       `yNNN     `NN");
            System.out.println( "dNo     /NNNNNm.       .dy.        .smNNNNNNNNs       NNNNNm-       `dNd     -Nm");
            System.out.println( "yNh     .NNNNNNd`       -Nmo.        -yNNNNNNNy       mNNNNNm.       .mo     +Ny");
            System.out.println( "/NN.     yNNNNNNy`       /NNm+`        :hNNNNNh       hNNNNNNh`       :.     hN/");
            System.out.println( "`mNs     -NNNNNNNo        oNNNd:`        /dNNNm       yNNNNNNNy             -Nm`");
            System.out.println( " +NN.     oNNNNNNN+        sNNNNy-        `/dNN       oNNNNNNNNo            hN+ ");

            System.out.println("---------------------------------------------------------------------------------------");
            System.out.println(model.getArtifactId());
            System.out.println("---------------------------------------------------------------------------------------");
            System.out.println("Version "+model.getVersion());
            System.out.println("From "+model.getGroupId());
            System.out.println("This part of "+model.getArtifactId()+" allows to compute distances between each pair fingerprint metabolites.");

            System.out.println();

            System.out.println("Required parameters");
            System.out.println("    -fingerprint <filepath> : a file with metabolite dbIdentifier");
            System.out.println("    -network <filepath> : a sbml file (metabolic network standard)");
            System.out.println("    -atommapping <filepath> : a file with Atom Atom Mapping");

            System.out.println();

            System.out.println("Facultative parameters");
            System.out.println("    -h : display help");
            System.out.println("    -matrixresult <filepath> : an output redirection of distance matrix");
            System.out.println("    -reactionresult <filepath> : an output redirection to keep path between metabolites");
            System.out.println("    -algo <filepath> : choose between ShortestAsUndirected and ValidShortest algorithm, default:ValidShortest");
            System.out.println("    -useMetExploreOutput : to return json for MetExplore integration");
            System.out.println("---------------------------------------------------------------------------------------");

            System.out.println("     ___  ___   ____   _____   ____  __    __  _____   _      _____   _____    ____  ");
            System.out.println("    /   |/   | | ___| |_   _| | ___| \\ \\  / / |  _  \\ | |    /  _  \\ |  _  \\  | ___| ");
            System.out.println("   / /|   /| | | |_     | |   | |_    \\ \\/ /  | |_| | | |    | | | | | |_| |  | |_   ");
            System.out.println("  / / |__/ | | |  _|    | |   |  _|    }  {   |  ___/ | |    | | | | |  _  /  |  _|  ");
            System.out.println(" / /       | | | |__    | |   | |__   / /\\ \\  | |     | |__  | |_| | | | \\ \\  | |__  ");
            System.out.println("/_/        |_| |____|   |_|   |____| /_/  \\_\\ |_|     |____| \\_____/ |_|  \\_\\ |____| ");
        }
        else
        {
            if(!this.useMetExploreOutput){
                //create a temp file
                File temp = File.createTempFile("closeErr", ".tmp");

                OutputStream output = new FileOutputStream(temp.getAbsolutePath());
                PrintStream printErr = new PrintStream(output);

                System.setErr(printErr);
            }
            String json = "";
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

                json = this.launchDistanceMatrix().toString();
            }
            else
            {
                json = this.launchDistanceMatrixFile().toString();
            }

            if(this.useMetExploreOutput){
                Gson gson = new GsonBuilder().setPrettyPrinting().create();
                JsonParser jp = new JsonParser();
                JsonElement je = jp.parse(json);
                String prettyJsonString = gson.toJson(je);

                System.out.println(prettyJsonString );
            }
        }
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
