import org.jdom2.Document;

/**
 * Created with IntelliJ IDEA.
 * Date: 03.01.14
 */
public abstract class Aligner extends Constants{

  int match;
  int mismatch;
  int gop;//gap open penalty
  int gep; //gap extension penalty

  public abstract Alignment align(byte[] s1, byte[] s2);

  public abstract int findMaxScore(byte[] s1, byte[] s2);

  public Aligner(Document document){}

}
