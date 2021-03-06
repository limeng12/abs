import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Map.Entry;

import htsjdk.samtools.*;


public class GetExon {

	public static void main(String[] args) throws Exception {
		
		String out_path=args[0];
		
		String bam_path=args[1];
		
		File bamFile = new File(bam_path);
		SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(bamFile);
		SAMRecordIterator r = sr.iterator();
		
		//File output_one= new File("output");
		BufferedWriter output_one = new BufferedWriter(new FileWriter(out_path) );

		HashMap<String,Integer>  all_exon= new HashMap<String,Integer>();
		
		while (r.hasNext()) {
			SAMRecord recode=r.next();
			//String att=recode.getStringAttribute("jI");
			
			int[] att;
			
			if(recode.hasAttribute("jI")) {
				att=recode.getSignedIntArrayAttribute("jI");
			}else {
				att=Reada5a3.get_splice_junction_pos(recode.getCigar(),recode.getStart() );
			}
			
			boolean is_neg=recode.getReadNegativeStrandFlag();
			
			if( (att.length<=2) || (att[0]==-1 )) {
				continue;
			}
			
			//for(int i=0;i<((att.length/2)-1);i++) {
				for(int i=0;i<(Math.floor(att.length/2)-1);i++) {

				String one_region=recode.getContig()+":"+(att[1+i*2]+1)+"-"+(att[2+i*2]-1);
				
				if(!all_exon.containsKey(one_region)) {
					all_exon.put(one_region, 1);
					continue;
				}
				
				int one_rcount=all_exon.get(one_region);
				all_exon.put(one_region, one_rcount+1);
				
			}
			
			
			/*	if(att.length>4) {
				String one_region=recode.getContig()+":"+(att[1]+1)+"-"+(att[2]-1);
				all_exon.add(one_region);
				
				
				one_region=recode.getContig()+":"+(att[3]+1)+"-"+(att[4]-1);
				all_exon.add(one_region);
				//System.out.println(  one_region );
				continue;
			}
			
			
			if(att.length>2) {
				String one_region=recode.getContig()+":"+(att[1]+1)+"-"+(att[2]-1);
				all_exon.add(one_region);
				//System.out.println(  one_region );
				continue;
			}*/
			

			
			
			
		}
		
        for(Entry<String, Integer> m :all_exon.entrySet()){
        	output_one.write(m.getKey()+"\t"+m.getValue());
			output_one.newLine();

        }
        output_one.close();
		
/*		it = all_exon.iterator(); // why capital "M"?
		while(it.hasNext()) {
			output_one.write(it.next());
			output_one.newLine();
		}
		output_one.close();*/
		
		
		r.close();
		sr.close();
		
	}

}
