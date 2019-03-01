import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

/**
 * Created by Gennady on 15.10.2016.
 */
public class FastaReader{

  protected static HashMap<String, String> read(String path, boolean replaceSpaces) throws IOException{
    HashMap<String, String> res = new HashMap<>();
    BufferedReader reader = new BufferedReader(new FileReader(path));
    StringBuilder builder = new StringBuilder();
    String name = null;
    for(String line = reader.readLine(); line != null; line = reader.readLine()){
      if(line.startsWith(">")){
        if(name != null){
          res.put(name, builder.toString());
        }
        name = line.substring(1);
        if(replaceSpaces){
          name = name.replace(' ', '_');
        }
        builder = new StringBuilder();
      }else{
        builder.append(line);
      }
    }
    if(name != null){
      res.put(name, builder.toString());
    }
    return res;
  }

  protected static HashMap<String, String> read(String path) throws IOException{
    return read(path, false);
  }

}
