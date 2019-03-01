import net.sourceforge.argparse4j.inf.Namespace;
import org.jdom2.Document;
import org.jdom2.Element;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

/**
 * Created by Геннадий on 06.12.2014.
 */
class DataReader extends Constants {
	
  private static final Pattern READ_NAME_TAIL = Pattern.compile("[/\\s].*");

  private String pathToReadsP1;
  private String pathToReadsP2;

  private String[] pathToReads1;
  private String[] pathToReads2;

  private static DataReader instance;

  public ArrayList<PairedRead> pairedReads;
  public int minReadLen;
  public int maxReadLen;

//  private DataReader(String pathToReadsP1, String pathToReadsP2){
//    this.pathToReadsP1 = pathToReadsP1;
//    this.pathToReadsP2 = pathToReadsP2;
//  }

  private void readData() throws IOException{
    if(pathToReads1 != null){
      pairedReads = readFilesWithReads(pathToReads1, pathToReads2);
    }else{
      pairedReads = readFilesWithReads(pathToReadsP1, pathToReadsP2);
    }
    minReadLen = Integer.MAX_VALUE;
    maxReadLen = 0;
    for(PairedRead read: pairedReads){
      if(!Arrays.equals(read.seq1, NULL_SEQ)){
        if(read.seq1.length < minReadLen){
          minReadLen = read.seq1.length;
        }
        if(read.seq1.length > maxReadLen){
          maxReadLen = read.seq1.length;
        }
      }
      if(!Arrays.equals(read.seq2, NULL_SEQ)){
        if(read.seq2.length < minReadLen){
          minReadLen = read.seq2.length;
        }
        if(read.seq2.length > maxReadLen){
          maxReadLen = read.seq2.length;
        }
      }
    }
  }

  private DataReader(Document document) throws IOException{
    Element element = document.getRootElement().getChild("Data");
    pathToReadsP1 = element.getChildText("pathToReads1");
    pathToReadsP2 = element.getChildText("pathToReads2");
    String[] paths1 = pathToReadsP1.split(",");
    String[] paths2 = pathToReadsP2.split(",");
    if(paths1.length > 1 || paths2.length > 1){
      pathToReads1 = paths1;
      pathToReads2 = paths1;
    }
    readData();
  }

  private DataReader(Document document, Namespace parsedArgs) throws Exception{
    pathToReadsP1 = parsedArgs.getString("pair1");
    pathToReadsP2 = parsedArgs.getString("pair2");
    if(pathToReadsP1 == null && pathToReadsP2 == null){
      Element element = document.getRootElement().getChild("Data");
      pathToReadsP1 = element.getChildText("pathToReads1");
      pathToReadsP2 = element.getChildText("pathToReads2");
    }else{
      if(pathToReadsP1 == null || pathToReadsP2 == null){
        throw new Exception("Both pair1 and pair2 should be specified or none of them!");
      }
    }
    String[] paths1 = pathToReadsP1.split(",");
    String[] paths2 = pathToReadsP2.split(",");
    if(paths1.length > 1 || paths2.length > 1){
      pathToReads1 = paths1;
      pathToReads2 = paths1;
    }
    readData();
  }

  public static DataReader getInstance(Document document){
    if(instance == null){
      try{
        instance = new DataReader(document);
      }catch(IOException e){
        e.printStackTrace();
      }
    }
    return instance;
  }

  public static DataReader getInstance(Document document, Namespace parsedArgs){
    if(instance == null){
      try{
        instance = new DataReader(document, parsedArgs);
      }catch(Exception e){
        e.printStackTrace();
      }
    }
    return instance;
  }

  private static BufferedReader newReadsReader(String path) throws IOException {
    InputStream is = new FileInputStream(path);
    if (path.toLowerCase().endsWith(".gz")) {
      is = new GZIPInputStream(is, 65536);
    }
    return new BufferedReader(new InputStreamReader(is), 32768);
  }

  private ArrayList<PairedRead> readFilesWithReads(String pathToReadsP1,
                                                     String pathToReadsP2)
      throws IOException {
    HashMap<String, PairedRead> readsMap = new HashMap<>();
    if (!pathToReadsP1.equals("-") && !pathToReadsP2.equals("-")) {
      BufferedReader readerP1 = newReadsReader(pathToReadsP1);
      String name;
      while ((name = readerP1.readLine()) != null){
        Matcher matcher = READ_NAME_TAIL.matcher(name);
        String readName;
        if(matcher.find()){
          readName = name.substring(1, matcher.start());
        }else{
          throw new IOException("File " + pathToReadsP1 + " have incorrect sequence identifier string");
        }
        byte[] sArr = readerP1.readLine().getBytes();
        PairedRead read = new PairedRead(readName, sArr, NULL_SEQ);
        read.q2 = NULL_SEQ;
        String qName = readerP1.readLine();
        if(!qName.startsWith("+") || (qName.length() > 1 && !qName.substring(1).equals(name.substring(1)))){
          throw new IOException("File " + pathToReadsP1 + " have incorrect quality string");
        }
        String q1 = readerP1.readLine();
        read.q1 = q1.getBytes();
        readsMap.put(readName, read);
      }
      BufferedReader readerP2 = newReadsReader(pathToReadsP2);
      while ((name = readerP2.readLine()) != null) {
        Matcher matcher = READ_NAME_TAIL.matcher(name);
        String readName;
        if(matcher.find()){
          readName = name.substring(1, matcher.start());
        }else{
          throw new IOException("File " + pathToReadsP2 + " have incorrect sequence identifier string");
        }
        PairedRead read = readsMap.get(readName);
        byte[] sArr = readerP2.readLine().getBytes();
        if (read == null) {
          read = new PairedRead(readName, NULL_SEQ, sArr);
          read.q1 = NULL_SEQ;
          readsMap.put(readName, read);
        } else {
          read.seq2 = sArr;
        }
        String qName = readerP2.readLine();
        if (!qName.startsWith("+") || (qName.length() > 1 && !qName.substring(1).equals(name.substring(1)))) {//
          throw new IOException("File " + pathToReadsP2 + " have incorrect quality string");
        }
        String q2 = readerP2.readLine();
        read.q2 = q2.getBytes();
      }
    }
    ArrayList<PairedRead> reads = new ArrayList<>(readsMap.size());
    reads.addAll(readsMap.values());
    return reads;
  }

  private ArrayList<PairedRead> readFilesWithReads(String[] pathToReads1,
                                                     String[] pathToReads2)
      throws IOException{
    HashMap<String, PairedRead> readsMap = new HashMap<>();
    for(String path: pathToReads1){
      BufferedReader readerP1 = newReadsReader(path);
      String name;
      while((name = readerP1.readLine()) != null){
        Matcher matcher = READ_NAME_TAIL.matcher(name);
        String readName;
        if(matcher.find()){
          readName = name.substring(1, matcher.start());
        }else{
          throw new IOException("File " + path + " have incorrect sequence identifier string");
        }
        byte[] sArr = readerP1.readLine().getBytes();
        PairedRead read = new PairedRead(readName, sArr, NULL_SEQ);
        read.q2 = NULL_SEQ;
        String qName = readerP1.readLine();
        if(!qName.startsWith("+") || (qName.length() > 1 && !qName.substring(1).equals(name.substring(1)))){
          throw new IOException("File " + path + " have incorrect quality string");
        }
        String q1 = readerP1.readLine();
        read.q1 = q1.getBytes();
        readsMap.put(readName, read);
      }
    }
    for(String path: pathToReads2){
      BufferedReader readerP2 = newReadsReader(path);
      String name;
      while ((name = readerP2.readLine()) != null) {
        Matcher matcher = READ_NAME_TAIL.matcher(name);
        String readName;
        if(matcher.find()){
          readName = name.substring(1, matcher.start());
        }else{
          throw new IOException("File " + path + " have incorrect sequence identifier string");
        }
        PairedRead read = readsMap.get(readName);
        byte[] sArr = readerP2.readLine().getBytes();
        if (read == null) {
          read = new PairedRead(readName, NULL_SEQ, sArr);
          read.q1 = NULL_SEQ;
          readsMap.put(readName, read);
        } else {
          read.seq2 = sArr;
        }
        String qName = readerP2.readLine();
        if (!qName.startsWith("+") || (qName.length() > 1 && !qName.substring(1).equals(name.substring(1)))) {//
          throw new IOException("File " + path + " have incorrect quality string");
        }
        String q2 = readerP2.readLine();
        read.q2 = q2.getBytes();
      }
    }
    ArrayList<PairedRead> reads = new ArrayList<>(readsMap.size());
    reads.addAll(readsMap.values());
    return reads;
  }



}
