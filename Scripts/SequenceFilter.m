% Reading the raw seq files 
Seq_delta = fastaread("Delta.fasta");
tb = struct2table(Seq_delta);
tb.len = zeros(length(tb.Sequence),1);
% Finding the length of each sequence
for i=1:height(tb)
    tb.len(i) = length(tb.Sequence{i});
end
% Plotting a histogram to see the distribution for sequence length
histogram(tb.len)
% Choosing the most common sequence length 
m = mode(tb.len);
% Filtering out these sequences 
TB = tb(tb.len == m,:);
% Comparing them with WT using BLAST for indels 
endseq = 'CCATGTGATTTTAATAGCTTCTTAGGAGAATGACAAAAAAAAAA';
endseq = length(endseq);
begseq = 'TTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTA';
begseq = length(begseq);
mutation_calibration = endseq+begseq;
TB_ready =table(TB.Header, TB.Sequence,'VariableNames',{'Header','Sequence'});
FilteredSeq = table2struct(TB_ready);
fastawrite("Delta-processed.fasta",FilteredSeq);

Ref_seq = TB.Sequence{1};
TB.nmut = zeros(height(TB),1);
for i=2:height(TB)
    TB.nmut(i) = length(TB.Sequence{i})-sum(TB.Sequence{i} == Ref_seq);
end