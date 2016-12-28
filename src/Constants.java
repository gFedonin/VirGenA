/**
 * Created with IntelliJ IDEA.
 * Date: 27.08.13
 */
public class Constants{

  static final byte GAP = (byte)'-';
  static final byte[] NULL_SEQ = {(byte)'*'};
  static final String NULL_SEQ_STR = "*";

  private static byte[] complement;
  static{
    complement = new byte[127];
    complement[(byte)'A'] = (byte)'T';
    complement[(byte)'T'] = (byte)'A';
    complement[(byte)'G'] = (byte)'C';
    complement[(byte)'C'] = (byte)'G';
    complement[(byte)GAP] = (byte)GAP;
    complement[(byte)'N'] = (byte)'N';
  }

  static byte[] lower;
  static{
    lower = new byte[127];
    lower[(byte)'A'] = (byte)'a';
    lower[(byte)'T'] = (byte)'t';
    lower[(byte)'G'] = (byte)'g';
    lower[(byte)'C'] = (byte)'c';
    lower[(byte)'N'] = (byte)'n';
  }

  static byte[] upper;
  static{
    upper = new byte[127];
    upper[(byte)'a'] = (byte)'A';
    upper[(byte)'t'] = (byte)'T';
    upper[(byte)'g'] = (byte)'G';
    upper[(byte)'c'] = (byte)'C';
    upper[(byte)'n'] = (byte)'N';
  }

  static byte[] getComplement(byte[] seq){
    int len = seq.length;
    byte[] res = new byte[len];
    for(int i = 0; i < len; i++){
      res[i] = complement[seq[len - i - 1]];
    }
    return res;
  }

  static int[] nToI;
  static byte[] iToN;
  static{
    nToI = new int['T' + 1];
    nToI[(byte)'A'] = 0;
    nToI[(byte)'T'] = 1;
    nToI[(byte)'G'] = 2;
    nToI[(byte)'C'] = 3;
    nToI[GAP] = 4;
    iToN = new byte[5];
    iToN[0] = (byte)'A';
    iToN[1] = (byte)'T';
    iToN[2] = (byte)'G';
    iToN[3] = (byte)'C';
    iToN[4] = GAP;
  }

  int anchorSize;

}
