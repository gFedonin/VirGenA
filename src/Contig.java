/**
 * Created with IntelliJ IDEA.
 * Date: 08.08.15
 */
class Contig{

  int alnStart;
  int alnEnd;
  int centroidStart;
  int centroidEnd;
  byte[] seq;
  byte[] seqAln;
  int[] coverage;
  int[] coverageAln;

  Contig(){}

  Contig(Contig contig){
    alnStart = contig.alnStart;
    alnEnd = contig.alnEnd;
    centroidStart = contig.centroidStart;
    centroidEnd = contig.centroidEnd;
    seq = contig.seq;
    seqAln = contig.seqAln;
    coverage = contig.coverage;
    coverageAln = contig.coverageAln;
  }

}
