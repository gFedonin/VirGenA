import htsjdk.samtools.*;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

/**
 * Created by Геннадий on 06.12.2014.
 */
class BamPrinter extends Constants{

  private int[] contigToRefIndex;


  private void setRecordFieldsMapped(MappedRead read, byte[] q, SAMRecord record,
                                     MappedRead mate, boolean isFirst, Reference genome){

    record.setReadString(new String(read.seq));
    int s = read.start + read.aln.start2;
    int refIndex = 0;
    int prev = 0;
    for(; refIndex < genome.contigEnds.length; refIndex++){
      int end = genome.contigEnds[refIndex];
      if(s < end){
        s -= prev;
        break;
      }
      prev = end;
    }
    record.setAlignmentStart(s + 1);
    setCigar(read, record);
    record.setMappingQuality(255);
    if(genome.isFragmented){
      record.setReferenceIndex(contigToRefIndex[refIndex]);
    }else{
      record.setReferenceIndex(refIndex);
    }
    record.setReadNegativeStrandFlag(read.reverse == 1);
    record.setReadUnmappedFlag(false);
    record.setReadPairedFlag(true);
    record.setFirstOfPairFlag(isFirst);
    record.setSecondOfPairFlag(!isFirst);
    if(read.reverse == 0){
      record.setBaseQualityString(new String(q));
    }else{
      record.setBaseQualityString(new String(invert(q)));
    }
    if(mate != null){
      record.setMateUnmappedFlag(false);
      int mateS = mate.start + mate.aln.start2;
      refIndex = 0;
      prev = 0;
      for(; refIndex < genome.contigEnds.length; refIndex++){
        int end = genome.contigEnds[refIndex];
        if(mateS < end){
          mateS -= prev;
          break;
        }
        prev = end;
      }
      record.setMateAlignmentStart(mateS + 1);
      record.setMateNegativeStrandFlag(mate.reverse == 1);
      if(genome.isFragmented){
        record.setMateReferenceIndex(contigToRefIndex[refIndex]);
      }else{
        record.setMateReferenceIndex(refIndex);
      }
    }else{
      record.setMateUnmappedFlag(true);
      record.setMateAlignmentStart(0);
      record.setMateNegativeStrandFlag(false);
      record.setMateReferenceName(NULL_SEQ_STR);
    }

  }



  private void setRecordFieldsUnmapped(byte[] seq, byte[] q, SAMRecord record,
                                       MappedRead mate, boolean isFirst, Reference genome){
    record.setCigarString(NULL_SEQ_STR);
    record.setAlignmentStart(0);
    record.setFirstOfPairFlag(isFirst);
    record.setSecondOfPairFlag(!isFirst);
    record.setReadString(new String(seq));
    record.setMappingQuality(0);
    record.setReferenceName(NULL_SEQ_STR);
    record.setReadPairedFlag(true);
    record.setReadUnmappedFlag(true);
    record.setBaseQualityString(new String(q));
    if(mate != null){
      record.setMateUnmappedFlag(false);
      int mateS = mate.start + mate.aln.start2;
      int refIndex = 0;
      int prev = 0;
      for(; refIndex < genome.contigEnds.length; refIndex++){
        int end = genome.contigEnds[refIndex];
        if(mateS < end){
          mateS -= prev;
          break;
        }
        prev = end;
      }
      record.setMateAlignmentStart(mateS + 1);
      record.setMateNegativeStrandFlag(mate.reverse == 1);
      if(genome.isFragmented){
        record.setMateReferenceIndex(contigToRefIndex[refIndex]);
      }else{
        record.setMateReferenceIndex(refIndex);
      }
    }else{
      record.setMateUnmappedFlag(true);
      record.setMateAlignmentStart(0);
      record.setMateNegativeStrandFlag(false);
      record.setMateReferenceName(NULL_SEQ_STR);
    }
  }


  private void setCigar(MappedRead read, SAMRecord record){
    Cigar cigar = new Cigar();
    if(read.aln.start1 > 0){
      CigarElement se = new CigarElement(read.aln.start1, CigarOperator.S);
      cigar.add(se);
    }
    CigarOperator operator = CigarOperator.M;
    int len = 0;
    int mNum = 0;
    for(int i = 0; i < read.aln.length; i++){
      byte c = read.aln.sequence1[i];
      if(c == GAP){
        if(operator == CigarOperator.D){
          len++;
        }else{
          CigarElement element = new CigarElement(len, operator);
          cigar.add(element);
          operator = CigarOperator.D;
          len = 1;
        }
      }else if(c < 'a'){
        if(operator == CigarOperator.M){
          len++;
        }else{
          CigarElement element = new CigarElement(len, operator);
          cigar.add(element);
          operator = CigarOperator.M;
          len = 1;
        }
        mNum++;
      }else{
        if(operator == CigarOperator.I){
          len++;
        }else{
          CigarElement element = new CigarElement(len, operator);
          cigar.add(element);
          operator = CigarOperator.I;
          len = 1;
        }
      }
    }
    CigarElement element = new CigarElement(len, operator);
    cigar.add(element);
    if(read.aln.end1 < read.seq.length){
      CigarElement se =
          new CigarElement(read.seq.length - read.aln.end1, CigarOperator.S);
      cigar.add(se);
    }
    record.setAttribute(ReservedTagConstants.XN, mNum - read.aln.identity);
    record.setCigar(cigar);
  }

  private byte[] invert(byte[] s){
    byte[] cArr = new byte[s.length];
    for(int i = 0; i < cArr.length; i++){
      cArr[i] = s[cArr.length - i - 1];
    }
    return cArr;
  }

  private static class SamRecordComparator implements Comparator<SAMRecord>{

    @Override
    public int compare(SAMRecord o1, SAMRecord o2){
      if(o1.getReferenceIndex().equals(o2.getReferenceIndex())){
        return o1.getAlignmentStart() - o2.getAlignmentStart();
      }else{
        return o1.getReferenceIndex() - o2.getReferenceIndex();
      }
    }

  }

  void printBAM(MappedData mappedData, Reference genome, String bamPath){
    if(bamPath.equals("-") || genome.length == 0){
      return;
    }
    SAMFileHeader header = new SAMFileHeader();
    SAMFileWriterFactory factory = new SAMFileWriterFactory();
    factory.setCreateIndex(true);
    header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
    contigToRefIndex = new int[genome.contigEnds.length];
    if(genome.isFragmented){
      int j = 0;
      int refID = 0;
      while(genome.fragmentEnds[j] == 0){
//        SAMSequenceRecord samSequenceRecord =
//            new SAMSequenceRecord(genome.fragmentNames[j], 0);
//        header.addSequence(samSequenceRecord);
//        refID ++;
        j ++;
      }
      for(int i = 0, k = 1, prev = 0; i < genome.contigEnds.length; i++){
        int end = genome.contigEnds[i];
        if(end == prev){
          continue;
        }
        SAMSequenceRecord samSequenceRecord;
        if(end == genome.fragmentEnds[j] && k == 1){
          samSequenceRecord =
              new SAMSequenceRecord(genome.fragmentNames[j], end - prev);
        }else{
          samSequenceRecord =
              new SAMSequenceRecord(genome.fragmentNames[j] + "_" + Integer.toString(k), end - prev);
          k++;
        }
        prev = end;
        header.addSequence(samSequenceRecord);
        contigToRefIndex[i] = refID;
        refID ++;
        if(end == genome.fragmentEnds[j]){
          j++;
          k = 1;
          while(j < genome.fragmentEnds.length - 1 && end == genome.fragmentEnds[j]){
//            samSequenceRecord =
//                new SAMSequenceRecord(genome.fragmentNames[j].replace(' ', '_'), 0);
//            header.addSequence(samSequenceRecord);
//            refID ++;
            j ++;
          }
        }
      }
    }else{
      if(genome.contigEnds.length == 1){
        SAMSequenceRecord samSequenceRecord = new SAMSequenceRecord(genome.name, genome.length);
        header.addSequence(samSequenceRecord);
      }else{
        for(int i = 0, prevCoord = 0; i < genome.contigEnds.length; i++){
          SAMSequenceRecord samSequenceRecord = new SAMSequenceRecord(
              genome.name + "_" + Integer.toString(i + 1), genome.contigEnds[i] - prevCoord);
          header.addSequence(samSequenceRecord);
          prevCoord = genome.contigEnds[i];
        }
      }
    }
    SAMFileWriter writer = factory.makeBAMWriter(header, true, new File(bamPath));
    ArrayList<SAMRecord> records = new ArrayList<>();
    ArrayList<SAMRecord> recordsUnaligned = new ArrayList<>();
    for(PairedRead pair : mappedData.concordant){
      SAMRecord record = new SAMRecord(header);
      record.setReadName(pair.name);
      record.setProperPairFlag(true);
      setRecordFieldsMapped(pair.r1, pair.q1, record, pair.r2, true, genome);
      records.add(record);
      record = new SAMRecord(header);
      record.setReadName(pair.name);
      record.setProperPairFlag(true);
      setRecordFieldsMapped(pair.r2, pair.q2, record, pair.r1, false, genome);
      records.add(record);
    }
    for(PairedRead pair : mappedData.discordant){
      SAMRecord record = new SAMRecord(header);
      record.setReadName(pair.name);
      record.setProperPairFlag(false);
      setRecordFieldsMapped(pair.r1, pair.q1, record, pair.r2, true, genome);
      records.add(record);
      record = new SAMRecord(header);
      record.setReadName(pair.name);
      record.setProperPairFlag(false);
      setRecordFieldsMapped(pair.r2, pair.q2, record, pair.r1, false, genome);
      records.add(record);
    }
    for(PairedRead pair : mappedData.leftMateMapped){
      SAMRecord record = new SAMRecord(header);
      record.setReadName(pair.name);
      record.setProperPairFlag(false);
      setRecordFieldsMapped(pair.r1, pair.q1, record, pair.r2, true, genome);
      records.add(record);
      record = new SAMRecord(header);
      record.setReadName(pair.name);
      record.setProperPairFlag(false);
      setRecordFieldsUnmapped(pair.seq2, pair.q2, record, pair.r1, false, genome);
      recordsUnaligned.add(record);
    }
    for(PairedRead pair : mappedData.rightMateMapped){
      SAMRecord record = new SAMRecord(header);
      record.setReadName(pair.name);
      record.setProperPairFlag(false);
      setRecordFieldsUnmapped(pair.seq1, pair.q1, record, pair.r2, true, genome);
      recordsUnaligned.add(record);
      record = new SAMRecord(header);
      record.setReadName(pair.name);
      record.setProperPairFlag(false);
      setRecordFieldsMapped(pair.r2, pair.q2, record, pair.r1, false, genome);
      records.add(record);
    }
    for(PairedRead pair : mappedData.unmapped){
      SAMRecord record = new SAMRecord(header);
      record.setReadName(pair.name);
      record.setProperPairFlag(false);
      setRecordFieldsUnmapped(pair.seq1, pair.q1, record, pair.r2, true, genome);
      recordsUnaligned.add(record);
      record = new SAMRecord(header);
      record.setReadName(pair.name);
      record.setProperPairFlag(false);
      setRecordFieldsUnmapped(pair.seq2, pair.q2, record, pair.r1, false, genome);
      recordsUnaligned.add(record);
    }
    Collections.sort(records, new SamRecordComparator());
    for(SAMRecord record: records){
      writer.addAlignment(record);
    }
    for(SAMRecord record: recordsUnaligned){
      writer.addAlignment(record);
    }
    writer.close();
  }

}
