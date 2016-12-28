import org.jdom2.Document;
import org.jdom2.input.SAXBuilder;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;

/**
 * Created with IntelliJ IDEA.
 * Date: 06.10.13
 */
public class MarkovModel extends Constants{

  private int K;
  private HashMap<String, float[]> freqs;
  private String[] genomes;

  MarkovModel(String[] genomes, int order){
    K = order;
    this.genomes = genomes;
    freqs = new HashMap<>();
    HashMap<String, int[]> counts = new HashMap<>();
    for(String genome: genomes){
      outer:
      for(int j = 0; j < genome.length() - K; j++){
        String s = genome.substring(j, j + K);
        for(int i = 0; i < K; i++){
          char c = s.charAt(i);
          if(c != 'A' && c != 'T' && c != 'G' && c != 'C'){
            continue outer;
          }
        }
        int[] cnt = counts.get(s);
        switch(genome.charAt(j + K)){
          case 'A':
            if(cnt == null){
              cnt = new int[4];
              counts.put(s, cnt);
            }
            cnt[0]++;
            break;
          case 'T':
            if(cnt == null){
              cnt = new int[4];
              counts.put(s, cnt);
            }
            cnt[1]++;
            break;
          case 'G':
            if(cnt == null){
              cnt = new int[4];
              counts.put(s, cnt);
            }
            cnt[2]++;
            break;
          case 'C':
            if(cnt == null){
              cnt = new int[4];
              counts.put(s, cnt);
            }
            cnt[3]++;
            break;
          default:
            //System.out.println("Unknown nucleotide at pos " + (j + K - 1));
        }
      }
      for(Map.Entry<String, int[]> entry : counts.entrySet()){
        float[] fr = new float[4];
        String key = entry.getKey();
        freqs.put(key, fr);
        int[] cnt = entry.getValue();
        int total = 0;
        for(int count : cnt){
          total += count;
        }
        float cf = 1.0f/total;
        fr[0] = cnt[0]*cf;
        fr[1] = fr[0] + cnt[1]*cf;
        fr[2] = fr[1] + cnt[2]*cf;
        fr[3] = 1.0f;
      }
    }
  }

  private static int findMaxK(String[] genomes, int maxK, float maxAvSumSqr){
    for(int k = 1; k <= maxK; k++){
      MarkovModel model = new MarkovModel(genomes, k);
      float sumSqr = 0;
      for(Map.Entry<String, float[]> entry: model.freqs.entrySet()){
        float[] freqs = entry.getValue();
        float prev = 0;
        for(float freq : freqs){
          float prob = freq - prev;
          sumSqr += prob*prob;
          prev = freq;
        }
      }
      if(sumSqr/model.freqs.size() > maxAvSumSqr){
        return k - 1;
      }
    }
    return maxK;
  }

  byte[] generate(int len){
    Random rnd = new Random();
    StringBuilder res = new StringBuilder();
    int genomeID = rnd.nextInt(genomes.length);
    String genome = genomes[genomeID];
    int genomePos = rnd.nextInt(genome.length() - K + 1);
    res.append(genome.substring(genomePos, genomePos + K));
    for(int j = K; j < len; j++){
      float[] fr = freqs.get(res.substring(j - K, j));
      if(fr == null){
        int n = rnd.nextInt(4);
        res.append(iToN[n]);
      }else{
        float pr = rnd.nextFloat();
        if(pr <= fr[0]){
          res.append('A');
        }else if(pr <= fr[1]){
          res.append('T');
        }else if(pr <= fr[2]){
          res.append('G');
        }else{
          res.append('C');
        }
      }
    }
    return res.toString().getBytes();
  }

  public static void main(String[] args){
    try{
      SAXBuilder jdomBuilder = new SAXBuilder();
      Document jdomDocument = jdomBuilder.build("./config.xml");
      ReferenceAlignment referenceAlignment = ReferenceAlignment.getInstance(jdomDocument);
      String[] genomes = new String[referenceAlignment.refAlns.size()];
      for(int i = 0; i < genomes.length; i++){
        genomes[i] = referenceAlignment.refAlns.get(i).seq;
      }
      System.out.println(findMaxK(genomes, 20, 0.65f));
      genomes = new String[1];
      genomes[0] = referenceAlignment.refAlns.get(0).seq;
      System.out.println(findMaxK(genomes, 10, 0.65f));
    }catch(Exception e){
      e.printStackTrace();
    }
  }

}
