import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by Gennady on 02.04.2016.
 */
class BlastParser{

  static ArrayList<BlastHit> parseBlastCSV(InputStream inputStream) throws IOException{
    ArrayList<BlastHit> res = new ArrayList<>();
    HashMap<String, BlastHit> map = new HashMap<>();
    BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
    for(String line = reader.readLine(); line != null; line = reader.readLine()){
      String[] fields = line.split(",");
      BlastHit hit = map.get(fields[1]);
      if(hit == null){
        hit = new BlastHit();
        hit.hitID = fields[1];
        map.put(fields[1], hit);
      }
      HSP hsp = new HSP();
      hsp.identity = Float.parseFloat(fields[2])/100;
      hsp.startQ = Integer.parseInt(fields[6]) - 1;
      hsp.endQ = Integer.parseInt(fields[7]);
      hsp.startS = Integer.parseInt(fields[8]) - 1;
      hsp.endS = Integer.parseInt(fields[9]);
      if(hsp.endS < hsp.startS){
        int end = hsp.startS + 1;
        hsp.startS = hsp.endS - 1;
        hsp.endS = end;
        hsp.reversed = true;
      }
      hsp.alnLen = Integer.parseInt(fields[3]);
      hsp.eValue = Float.parseFloat(fields[10]);
      hsp.score = Float.parseFloat(fields[11]);
      hit.hsps.add(hsp);
    }
    res.addAll(map.values());
    return res;
  }

}
