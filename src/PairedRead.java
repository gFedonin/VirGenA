/**
 * Created by Gennady
 */
class PairedRead{

  String name;
  byte[] seq1;
  byte[] seq2;
  byte[] q1;
  byte[] q2;
  MappedRead r1;
  MappedRead r2;

  PairedRead(String name, byte[] seq1, byte[] seq2){
    this.name = name;
    this.seq1 = seq1;
    this.seq2 = seq2;
  }

  PairedRead(String name, byte[] seq1, byte[] seq2, byte[] q1, byte[] q2){
    this.name = name;
    this.seq1 = seq1;
    this.seq2 = seq2;
    this.q1 = q1;
    this.q2 = q2;
  }

  PairedRead(Read read1, Read read2){
    name = read1.name;
    if(read1.n == '1'){
      seq1 = read1.seq;
      q1 = read1.q;
      seq2 = read2.seq;
      q2 = read2.q;
    }else{
      seq1 = read2.seq;
      q1 = read2.q;
      seq2 = read1.seq;
      q2 = read1.q;
    }
  }

  PairedRead(){

  }

  PairedRead(PairedRead pairedRead){
    this.name = pairedRead.name;
    this.seq1 = pairedRead.seq1;
    this.seq2 = pairedRead.seq2;
    this.q1 = pairedRead.q1;
    this.q2 = pairedRead.q2;
    this.r1 = pairedRead.r1;
    this.r2 = pairedRead.r2;
  }

}
