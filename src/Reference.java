import gnu.trove.list.array.TIntArrayList;
import org.jdom2.Document;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by Геннадий on 14.02.2015.
 */
public class Reference extends Constants{

  public String name;
  public String seq;
  public String aln;
  private int[] alnToSeq;
  int[] seqToAln;
  HashMap<String, int[]> index;
  ArrayList<PairedRead> reads;
  String consensus;
  int refID;
  int[] contigEnds;
  int[] fragmentEnds;
  boolean isFragmented;
  String[] fragmentNames;
  byte[] seqB;
  int length;

  public Reference(String name, String aln){
    this.name = name;
    this.aln = aln;
    alnToSeq = new int[aln.length()];
    int strLen = 0;
    for(int i = 0; i < alnToSeq.length; i++){
      if(aln.charAt(i) != '-'){
        alnToSeq[i] = strLen;
        strLen ++;
      }else{
        alnToSeq[i] = -1;
      }
    }
    seqToAln = new int[strLen];
    char[] s = new char[strLen];
    for(int i = 0, j = 0; i < alnToSeq.length; i++){
      if(alnToSeq[i] != -1){
        seqToAln[j] = i;
        s[j] = aln.charAt(i);
        j ++;
      }
    }
    seq = new String(s);
  }

  int alnToSeqStart(int start){
    int res = alnToSeq[start];
    if(res == -1){
      int s = start + 1;
      while(s < alnToSeq.length && alnToSeq[s] == -1){
        s++;
      }
      if(s == alnToSeq.length){
        return  -1;
      }
      return alnToSeq[s];
    }else{
      return res;
    }
  }

  int alnToSeqEnd(int end){
    int res = alnToSeq[end];
    if(res == -1){
      int e = end - 1;
      while(e >= 0 && alnToSeq[e] == -1){
        e--;
      }
      if(e == -1){
        return  -1;
      }
      return alnToSeq[e] + 1;
    }else{
      return res;
    }
  }

//  public Reference(String seq){
//    this.seqB = seq.getBytes();
//    length = seq.length();
//  }

  public Reference(Document document){
    String pathToGenome = document.getRootElement().getChildText("Reference");
    int threadNum = Integer.parseInt(document.getRootElement().getChildText("ThreadNumber"));
    try{
      readGenome(pathToGenome);
      index = StringIndexer.buildIndexPara(seq, KMerCounter.getInstance(document).K, threadNum);
    }catch(Exception e){
      e.printStackTrace();
    }
  }

  public void buildIndex(Document document){
    length = seq.length();
    seqB = seq.getBytes();
    int threadNum = Integer.parseInt(document.getRootElement().getChildText("ThreadNumber"));
    try{
      index = StringIndexer.buildIndexPara(seq, KMerCounter.getInstance(document).K, threadNum);
    }catch (Exception e){
      e.printStackTrace();
    }
  }

  public Reference(ConsensusBuilderWithReassembling consensusBuilder){
    name = consensusBuilder.genomeName;
    seq = consensusBuilder.finalConsensus;
    seqB = seq.getBytes();
    length = seq.length();
    isFragmented = consensusBuilder.genomeIsFragmented;
    fragmentNames = consensusBuilder.fragmentNames;
    fragmentEnds = consensusBuilder.fragmentEnds;
    contigEnds = consensusBuilder.contigEnds;
    try{
      index = StringIndexer.buildIndexPara(seq, consensusBuilder.K, consensusBuilder.threadNum);
    }catch (Exception e){
      e.printStackTrace();
    }
  }

  public Reference(String name, String g, int K, int[] contigEnds, int[] fragmentEnds,
                String[] fragmentNames, boolean isFragmented, int threadNum){
    this.name = name;
    seq = g;
    this.contigEnds = contigEnds;
    this.fragmentEnds = fragmentEnds;
    this.fragmentNames = fragmentNames;
    this.isFragmented = isFragmented;
    length = seq.length();
    seqB = seq.getBytes();
    try{
      index = StringIndexer.buildIndexPara(seq, K, threadNum);
    }catch (Exception e){
      e.printStackTrace();
    }
  }

  public Reference(byte[] g, int K, Reference genome, int threadNum){
    this(g, K, genome, true, threadNum);
  }

  public Reference(byte[] g, int K, Reference genome, boolean buildIndex, int threadNum){
    name = genome.name;
    seqB = g;
    seq = new String(seqB);
    length = seq.length();
    contigEnds = genome.contigEnds;
    fragmentEnds = genome.fragmentEnds;
    isFragmented = genome.isFragmented;
    fragmentNames = genome.fragmentNames;
    if(buildIndex){
      try{
        index = StringIndexer.buildIndexPara(seq, K, threadNum);
      }catch(Exception e){
        e.printStackTrace();
      }
    }
  }

  private void readGenome(String pathToGenome) throws IOException{
    BufferedReader reader = new BufferedReader(new FileReader(pathToGenome));
    String line;
    StringBuilder builder = new StringBuilder();
    TIntArrayList ends = new TIntArrayList();
    ArrayList<String> names = new ArrayList<>();
    int prevLen = 0;
    while((line = reader.readLine()) != null){
      if(line.startsWith(">")){
        names.add(line.substring(1));
        if(prevLen > 0){
          ends.add(prevLen);
        }
      }else{
        builder.append(line);
        prevLen += line.length();
      }
    }
    ends.add(prevLen);
    contigEnds = ends.toArray();
    fragmentEnds = contigEnds;
    seq = builder.toString();
    length = seq.length();
    seqB = seq.getBytes();
    fragmentNames = names.toArray(new String[names.size()]);
    if(fragmentNames.length > 1){
      isFragmented = true;
      int i = pathToGenome.lastIndexOf(".fasta");
      if(i == -1){
        i = pathToGenome.lastIndexOf(".fa");
      }
      int j = pathToGenome.lastIndexOf('/');
      if(j == -1){
        if(i == -1){
          name = pathToGenome;
        }else{
          name = pathToGenome.substring(0, i);
        }
      }else{
        if(i == -1){
          name = pathToGenome.substring(j + 1);
        }else{
          name = pathToGenome.substring(j + 1, i);
        }
      }
    }else{
      name = fragmentNames[0];
    }
  }

  public Reference(String pathToGenome, int K, int threadNum){
    try{
      readGenome(pathToGenome);
      index = StringIndexer.buildIndexPara(seq, K, threadNum);
    }catch(Exception e){
      e.printStackTrace();
    }
  }


}
