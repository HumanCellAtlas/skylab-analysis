# Manual Cellranger runs

The following lists the raw commands issued to perform perform count matrix generation with cellranger version 3.0.2. To reproduce this download cellranger 3.0.2 and amend the paths below.

```
PATH=/home/nbarkas/disk3/cellranger/cellranger-3.0.2:$PATH

cellranger count --id t_4k --fastqs /home/nbarkas/disk3/4k_panT/data --sample t_4k --transcriptome /home/nbarkas/disk3/cellranger/GRCh38 --localcores 24 --localmem 75

cellranger count --id pbmc4k --fastqs /home/nbarkas/disk3/4k_pbmc/data --sample pbmc4k --transcriptome /home/nbarkas/disk3/cellranger/GRCh38 --localcores 24 --localmem 75

cellranger count --id pbmc8k --fastqs /home/nbarkas/disk3/8k_pbmc/data --sample pbmc8k --transcriptome /home/nbarkas/disk3/cellranger/GRCh38 --localcores 24 --localmem 75

cellranger count --id HCATisStabAug177376561 --fastqs /home/nbarkas/disk3/ischaemia_1/data --sample HCATisStabAug177376561 --transcriptome /home/nbarkas/disk3/cellranger/GRCh38 --localcores 24 --localmem 75

cellranger count --id MantonBM1_HiSeq_1 --fastqs /home/nbarkas/disk3/MantonBM1_1/data --sample MantonBM1_HiSeq_1 --transcriptome /home/nbarkas/disk3/cellranger/GRCh38 --localcores 24 --localmem 75
```