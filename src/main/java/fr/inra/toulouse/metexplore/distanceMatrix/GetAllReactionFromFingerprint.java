package fr.inra.toulouse.metexplore.distanceMatrix;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonElement;
import com.google.gson.JsonParser;
import fr.inra.toulouse.metexplore.met4j_core.biodata.*;
import org.json.simple.JSONArray;
import org.json.simple.JSONObject;
import org.json.simple.parser.JSONParser;
import org.json.simple.parser.ParseException;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import javax.xml.parsers.ParserConfigurationException;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;

/**
 * The Class MeSHAnalysis.
 */
public class GetAllReactionFromFingerprint {

    private JSONObject lightestpath;
    /**
     * metabolic network
     */
    @Option(name = "-fingerprint", usage = "Global network", metaVar = "fingerprintPath", required = true)
    private String fingerprintPath = "fingerprintPath";

    /**
     * metabolic network
     */
    @Option(name = "-lightestpath", usage = "lightestpath", metaVar = "lightestpath", required = true)
    private String lightestpathFile = "lightestpath.json";

    /**
     * path length threshold
     *
     */
    @Option(name = "-threshold", usage = "threshold", metaVar = "threshold", required = false)
    private int threshold = 100;

    /**
     * Constructor
     */
    private GetAllReactionFromFingerprint() {
        super();
    }
    public GetAllReactionFromFingerprint(String fingerprintPath, String lightestpathFile) {
        super();
        this.fingerprintPath = fingerprintPath;
        this.lightestpathFile = lightestpathFile;
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
        GetAllReactionFromFingerprint app = new GetAllReactionFromFingerprint();

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


        System.err.println("Started");
        JSONParser parser = new JSONParser();
        JSONObject networkArg = null;
        try {
            this.lightestpath = (JSONObject) parser.parse(new FileReader(this.lightestpathFile));
        } catch (ParseException e) {
            e.printStackTrace();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return this.launch();
    }

    public JSONObject launch(){
        /**
         * Get metabolic profile as a tabular file with one column with metabolite identifiers used in the network
         */
        HashSet<String> metabolites=new HashSet<String>();
        try {
            BufferedReader reader=new BufferedReader(new FileReader(this.fingerprintPath));
            try {
                String line=reader.readLine();
                while(line!=null)
                {
                    String[] elements=line.split("\t");
                    metabolites.add(elements[0]);
                    line=reader.readLine();
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        } catch (FileNotFoundException e) {
            System.err.println("Please provide fingerprint file");
            e.printStackTrace();
        }
        System.err.println(metabolites.toString());
        HashSet<String> reactionsList = new HashSet<String>();

        for(String source : metabolites){
            for(String target : metabolites){

                int compare = source.compareTo(target);
                String pathName;
                if (compare < 0){
                    pathName = source+target;
                }
                else
                {
                    pathName = target+source;
                }

                if(this.lightestpath.keySet().contains(pathName)){
                    JSONArray reactions = (JSONArray) this.lightestpath.get(pathName);
                    if(reactions.size()<=threshold)
                    {
                        for (int i = 0; i < reactions.size(); ++i) {
                            String reac = (String) reactions.get(i);
                            reactionsList.add(reac);
                        }
                    }
                }
            }
        }
        for (String s : reactionsList) {
            System.err.println(s);
        }
        JSONArray json_results_reactionslist = new JSONArray();
        for (String rea : reactionsList) {
            json_results_reactionslist.add(rea);
        }


        JSONObject jsonObject = new JSONObject();
        jsonObject.put("message", new String("Results are displayed in file."));
        jsonObject.put("success", "true");
        jsonObject.put("reactions", json_results_reactionslist);


        System.err.println("Finished");

        return jsonObject;

    }
}
