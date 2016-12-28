/**
 * Created by Gennady on 02.04.2016.
 */
public class HSP implements Comparable{

  public float identity;
  public float score;
  public int startQ;
  public int endQ;
  public int startS;
  public int endS;
  public int alnLen;
  public float eValue;
  public boolean reversed;

  @Override
  public int compareTo(Object o){
    HSP hsp = (HSP) o;
    return startQ - hsp.startQ;
  }
}
