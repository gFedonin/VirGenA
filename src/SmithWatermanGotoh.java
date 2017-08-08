import org.jdom2.Document;
import org.jdom2.Element;
import org.jdom2.input.SAXBuilder;

class SmithWatermanGotoh extends Aligner{

    SmithWatermanGotoh(Document document){
        super(document);
        Element element = document.getRootElement().getChild("Mapper").getChild("Aligner");
        match = Integer.parseInt(element.getChildText("Match"));
        mismatch = Integer.parseInt(element.getChildText("Mismatch"));
        gop = Integer.parseInt(element.getChildText("GapOpenPenalty"));
        gep = Integer.parseInt(element.getChildText("GapExtensionPenalty"));
    }

    private class Cell{
        /**
         * Row of the cell
         */
        public short row;
        /**
         * Column of the cell
         */
        public short col;
        /**
         * Alignment score at this cell
         */
        public float score;
    }

    private class Directions{
        /**
         * Traceback direction stop
         */
        public static final byte STOP = 0;
        /**
         * Traceback direction left
         */
        public static final byte LEFT = 1;
        /**
         * Traceback direction diagonal
         */
        public static final byte DIAGONAL = 2;
        /**
         * Traceback direction up
         */
        public static final byte UP = 3;
    }

    @Override
    public Alignment align(byte[] s1, byte[] s2){

        int m = s1.length + 1;
        int n = s2.length + 1;

        byte[] pointers = new byte[m*n];

        // Initializes the boundaries of the traceback matrix to STOP.
        for(int i = 0, k = 0; i < m; i++, k += n){
            pointers[k] = Directions.STOP;
        }
        for(int j = 1; j < n; j++){
            pointers[j] = Directions.STOP;
        }

        short[] sizesOfVerticalGaps = new short[m*n];
        short[] sizesOfHorizontalGaps = new short[m*n];
        for(int i = 0, k = 0; i < m; i++, k += n){
            for(int j = 0; j < n; j++){
                sizesOfVerticalGaps[k + j] = sizesOfHorizontalGaps[k + j] = 1;
            }
        }

        Cell cell = construct(s1, s2, pointers, sizesOfVerticalGaps, sizesOfHorizontalGaps);
        return traceback(s1, s2, pointers, cell, sizesOfVerticalGaps, sizesOfHorizontalGaps);
    }


    private Cell construct(byte[] s1, byte[] s2, byte[] pointers, short[] sizesOfVerticalGaps,
                           short[] sizesOfHorizontalGaps){

        short m = (short) (s1.length + 1);
        short n = (short) (s2.length + 1);

        float f; // score of alignment x1...xi to y1...yi if xi aligns to yi
        float[] g = new float[n]; // score if xi aligns to a gap after yi
        float h; // score if yi aligns to a gap after xi
        float[] v = new float[n]; // best score of alignment x1...xi to y1...yi
        float vDiagonal;

        g[0] = Float.NEGATIVE_INFINITY;
        h = Float.NEGATIVE_INFINITY;
        v[0] = 0;

        for(int j = 1; j < n; j++){
            g[j] = Float.NEGATIVE_INFINITY;
            v[j] = 0;
        }

        float similarityScore, g1, g2, h1, h2;

        Cell cell = new Cell();

        for(int i = 1, k = n; i < m; i++, k += n){
            h = Float.NEGATIVE_INFINITY;
            vDiagonal = v[0];
            for(int j = 1, l = (k + 1); j < n; j++, l++){
                similarityScore = (s1[i - 1] == s2[j - 1]) ? match : mismatch;

                // Fill the matrices
                f = vDiagonal + similarityScore;

                g1 = g[j] - gep;
                g2 = v[j] - gop;
                if(g1 > g2){
                    g[j] = g1;
                    sizesOfVerticalGaps[l] = (short) (sizesOfVerticalGaps[l - n] + 1);
                }else{
                    g[j] = g2;
                }

                h1 = h - gep;
                h2 = v[j - 1] - gop;
                if(h1 > h2){
                    h = h1;
                    sizesOfHorizontalGaps[l] = (short) (sizesOfHorizontalGaps[l - 1] + 1);
                }else{
                    h = h2;
                }

                vDiagonal = v[j];
                v[j] = maximum(f, g[j], h, 0);

                // Determine the traceback direction
                if(v[j] == 0){
                    pointers[l] = Directions.STOP;
                }else if(v[j] == f){
                    pointers[l] = Directions.DIAGONAL;
                }else if(v[j] == g[j]){
                    pointers[l] = Directions.UP;
                }else{
                    pointers[l] = Directions.LEFT;
                }

                // Set the traceback start at the current cell i, j and identity
                if(v[j] > cell.score){
                    cell.row = (short) i;
                    cell.col = (short) j;
                    cell.score = v[j];
                }
            }
        }
        return cell;
    }


    private Alignment traceback(byte[] a1, byte[] a2,
                                byte[] pointers, Cell cell, short[] sizesOfVerticalGaps,
                                short[] sizesOfHorizontalGaps){


        int n = a2.length + 1;

        Alignment alignment = new Alignment();
        alignment.score = Math.round(cell.score);

        int maxlen = a1.length + a2.length; // maximum length after the
        // aligned sequences

        byte[] reversed1 = new byte[maxlen]; // reversed sequence #1

        int len1 = 0; // length of sequence #1 after alignment

        short identity = 0; // count of identitcal pairs

        byte c1, c2;

        int i = cell.row; // traceback start row
        alignment.end1 = i;
        int j = cell.col; // traceback start col
        alignment.end2 = j;
        int k = i*n;

        boolean stillGoing = true; // traceback flag: true -> continue & false
        // -> stop

        while(stillGoing){
            switch(pointers[k + j]){
                case Directions.UP:
                    for(int l = 0, len = sizesOfVerticalGaps[k + j]; l < len; l++){
                        reversed1[len1++] = lower[a1[--i]];
                        k -= n;
                    }
                    break;
                case Directions.DIAGONAL:
                    c1 = a1[--i];
                    c2 = a2[--j];
                    k -= n;
                    reversed1[len1++] = c1;
                    if(c1 == c2){
                        identity++;
                    }
                    break;
                case Directions.LEFT:
                    for(int l = 0, len = sizesOfHorizontalGaps[k + j]; l < len; l++){
                        reversed1[len1++] = GAP;
                        --j;
                    }
                    break;
                case Directions.STOP:
                    stillGoing = false;
            }
        }

        alignment.sequence1 = reverse(reversed1, len1);
        alignment.length = alignment.sequence1.length;
        alignment.start1 = i;
        alignment.start2 = j;
        alignment.identity = identity;

        return alignment;
    }


    private float maximum(float a, float b, float c, float d){
        if(a > b){
            if(a > c){
                return a > d ? a : d;
            }else{
                return c > d ? c : d;
            }
        }else if(b > c){
            return b > d ? b : d;
        }else{
            return c > d ? c : d;
        }
    }

    @Override
    public int findMaxScore(byte[] s1, byte[] s2){

        short m = (short) (s1.length + 1);
        short n = (short) (s2.length + 1);

        float f; // score of alignment x1...xi to y1...yi if xi aligns to yi
        float[] g = new float[n]; // score if xi aligns to a gap after yi
        float h; // score if yi aligns to a gap after xi
        float[] v = new float[n]; // best score of alignment x1...xi to y1...yi
        float vDiagonal;

        g[0] = Float.NEGATIVE_INFINITY;
        h = Float.NEGATIVE_INFINITY;
        v[0] = 0;

        for(int j = 1; j < n; j++){
            g[j] = Float.NEGATIVE_INFINITY;
            v[j] = 0;
        }

        float similarityScore, g1, g2, h1, h2;

        float maxScore = 0;

        for(short i = 1, k = n; i < m; i++, k += n){
            h = Float.NEGATIVE_INFINITY;
            vDiagonal = v[0];
            for(short j = 1, l = (short) (k + 1); j < n; j++, l++){
                similarityScore = (s1[i - 1] == s2[j - 1]) ? match : mismatch;

                // Fill the matrices
                f = vDiagonal + similarityScore;

                g1 = g[j] - gep;
                g2 = v[j] - gop;
                if(g1 > g2){
                    g[j] = g1;
                }else{
                    g[j] = g2;
                }

                h1 = h - gep;
                h2 = v[j - 1] - gop;
                if(h1 > h2){
                    h = h1;
                }else{
                    h = h2;
                }

                vDiagonal = v[j];
                v[j] = maximum(f, g[j], h, 0);


                // Set the traceback start at the current cell i, j and identity
                if(v[j] > maxScore){
                    maxScore = v[j];
                }
            }
        }
        return Math.round(maxScore);
    }

    /**
     * Reverses an array of chars
     *
     * @param a
     * @param len
     * @return the input array of char reserved
     */
    private byte[] reverse(byte[] a, int len){
        byte[] b = new byte[len];
        for(int i = len - 1, j = 0; i >= 0; i--, j++){
            b[j] = a[i];
        }
        return b;
    }

    public static void main(String[] args){
        try{
            SAXBuilder jdomBuilder = new SAXBuilder();
            Document jdomDocument = jdomBuilder.build("./config.xml");
            SmithWatermanGotoh aligner = new SmithWatermanGotoh(jdomDocument);
            String read = "AATATAGGAAAATAAGAAAACAAAAGAAAATAAACAGGTTAATTGATAGAATAAGAGAAAGAGCAGAAGACAGTGGCAATGAGAGCGAGGGAGATCAAGAAGAATTGTCAGCACTGGTGGTGGAGATGGGGCATCATGCTCCTTGGGATGTGGATGATTTGTAGTGCTGCAGAACAATTGTGGGTTACAGTTTATTATGGGGTTCCTGTGTGGAGAGATGCAGATACCACCCTATTTTGTGCATCAGATG";
            ReferenceAlignment refAlignment = ReferenceAlignment.getInstance(jdomDocument);
            Alignment aln = aligner.align(read.getBytes(), refAlignment.refSeqs.get("B.US.86.JRFL_JR_FL.U63632").seqB);
            System.out.printf("Score1 B = %d start = %d end = %d\n", aln.score, aln.start2, aln.end2);
            aln = aligner.align(read.getBytes(), refAlignment.refSeqs.get("01_AE.TH.90.CM240.U54771").seqB);
            System.out.printf("Score1 01_AE = %d start = %d end = %d\n", aln.score, aln.start2, aln.end2);
            byte[] readB = getComplement("GTGGCCCAGATATTGTGCATTTCTGTCTCATGTCCTTTAGCATCTGATGCACAAAATAGGGTGGTATCTGCATCTCTCCACACAGGAACCCCATAATAAACTGTAACCCACAATTGTTCTGCAGCACTACAAATCATCCACATCCCAAGGAGCATGATGCCCCATCTCCACCACCAGTGCTGACAATTCTTCTTGATCTCCCTCGCTCTCATTGCCACTGTCTTCTGCTCTTTCTCTTATTCTATCAATT".getBytes());
            aln = aligner.align(readB, refAlignment.refSeqs.get("B.US.86.JRFL_JR_FL.U63632").seqB);
            System.out.printf("Score2 B = %d start = %d end = %d\n", aln.score, aln.start2, aln.end2);
            aln = aligner.align(readB, refAlignment.refSeqs.get("01_AE.TH.90.CM240.U54771").seqB);
            System.out.printf("Score2 01_AE = %d start = %d end = %d\n", aln.score, aln.start2, aln.end2);
        }catch(Exception e){
            e.printStackTrace();
        }
    }

}