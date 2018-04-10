package fr.inra.toulouse.metexplore;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.ProtocolException;
import java.net.URL;
import java.util.HashSet;

public class MetaboWorkBenchAccess {

	public static void main(String[] args){
		System.out.println("start access");
		System.out.println(getFingerprint());
	}
	
	protected static HashSet<String> getFingerprint(){
		HashSet<String> fingerprint=new HashSet<String>();
		try {
			URL url = new URL("http://www.metabolomicsworkbench.org/rest/compound/regno/11/all");
			HttpURLConnection conn = (HttpURLConnection) url.openConnection();
			conn.setRequestMethod("GET");
			conn.setRequestProperty("Accept", "application/json");
			BufferedReader br = new BufferedReader(new InputStreamReader(
					(conn.getInputStream())));
			String output;
			System.out.println("Output from Server .... \n");
			while ((output = br.readLine()) != null) {
				System.out.println(output);
			}

			conn.disconnect();
			
		} catch (MalformedURLException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ProtocolException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		
		return fingerprint;
		
	}
}
