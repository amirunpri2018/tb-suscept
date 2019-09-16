#!/bin/bash
set -eu

# Download GSE94116 FASTQ files from SRA with sra-tools

echo "Downloading SRR5206790: r01-resist-infect (1 of 50)"
fastq-dump --stdout --gzip SRR5206790 > r01-resist-infect.fastq.gz
echo "Downloading SRR5206791: r01-resist-noninf (2 of 50)"
fastq-dump --stdout --gzip SRR5206791 > r01-resist-noninf.fastq.gz
echo "Downloading SRR5206792: r02-resist-infect (3 of 50)"
fastq-dump --stdout --gzip SRR5206792 > r02-resist-infect.fastq.gz
echo "Downloading SRR5206793: r02-resist-noninf (4 of 50)"
fastq-dump --stdout --gzip SRR5206793 > r02-resist-noninf.fastq.gz
echo "Downloading SRR5206794: r03-resist-infect (5 of 50)"
fastq-dump --stdout --gzip SRR5206794 > r03-resist-infect.fastq.gz
echo "Downloading SRR5206795: r03-resist-noninf (6 of 50)"
fastq-dump --stdout --gzip SRR5206795 > r03-resist-noninf.fastq.gz
echo "Downloading SRR5206796: r04-resist-infect (7 of 50)"
fastq-dump --stdout --gzip SRR5206796 > r04-resist-infect.fastq.gz
echo "Downloading SRR5206797: r04-resist-noninf (8 of 50)"
fastq-dump --stdout --gzip SRR5206797 > r04-resist-noninf.fastq.gz
echo "Downloading SRR5206798: r05-resist-infect (9 of 50)"
fastq-dump --stdout --gzip SRR5206798 > r05-resist-infect.fastq.gz
echo "Downloading SRR5206799: r05-resist-noninf (10 of 50)"
fastq-dump --stdout --gzip SRR5206799 > r05-resist-noninf.fastq.gz
echo "Downloading SRR5206800: r06-resist-infect (11 of 50)"
fastq-dump --stdout --gzip SRR5206800 > r06-resist-infect.fastq.gz
echo "Downloading SRR5206801: r06-resist-noninf (12 of 50)"
fastq-dump --stdout --gzip SRR5206801 > r06-resist-noninf.fastq.gz
echo "Downloading SRR5206802: r07-resist-infect (13 of 50)"
fastq-dump --stdout --gzip SRR5206802 > r07-resist-infect.fastq.gz
echo "Downloading SRR5206803: r07-resist-noninf (14 of 50)"
fastq-dump --stdout --gzip SRR5206803 > r07-resist-noninf.fastq.gz
echo "Downloading SRR5206804: r08-resist-infect (15 of 50)"
fastq-dump --stdout --gzip SRR5206804 > r08-resist-infect.fastq.gz
echo "Downloading SRR5206805: r08-resist-noninf (16 of 50)"
fastq-dump --stdout --gzip SRR5206805 > r08-resist-noninf.fastq.gz
echo "Downloading SRR5206806: r09-resist-infect (17 of 50)"
fastq-dump --stdout --gzip SRR5206806 > r09-resist-infect.fastq.gz
echo "Downloading SRR5206807: r09-resist-noninf (18 of 50)"
fastq-dump --stdout --gzip SRR5206807 > r09-resist-noninf.fastq.gz
echo "Downloading SRR5206808: r10-resist-infect (19 of 50)"
fastq-dump --stdout --gzip SRR5206808 > r10-resist-infect.fastq.gz
echo "Downloading SRR5206809: r10-resist-noninf (20 of 50)"
fastq-dump --stdout --gzip SRR5206809 > r10-resist-noninf.fastq.gz
echo "Downloading SRR5206810: r11-resist-infect (21 of 50)"
fastq-dump --stdout --gzip SRR5206810 > r11-resist-infect.fastq.gz
echo "Downloading SRR5206811: r11-resist-noninf (22 of 50)"
fastq-dump --stdout --gzip SRR5206811 > r11-resist-noninf.fastq.gz
echo "Downloading SRR5206812: r12-resist-infect (23 of 50)"
fastq-dump --stdout --gzip SRR5206812 > r12-resist-infect.fastq.gz
echo "Downloading SRR5206813: r12-resist-noninf (24 of 50)"
fastq-dump --stdout --gzip SRR5206813 > r12-resist-noninf.fastq.gz
echo "Downloading SRR5206814: r13-resist-infect (25 of 50)"
fastq-dump --stdout --gzip SRR5206814 > r13-resist-infect.fastq.gz
echo "Downloading SRR5206815: r13-resist-noninf (26 of 50)"
fastq-dump --stdout --gzip SRR5206815 > r13-resist-noninf.fastq.gz
echo "Downloading SRR5206816: r14-resist-infect (27 of 50)"
fastq-dump --stdout --gzip SRR5206816 > r14-resist-infect.fastq.gz
echo "Downloading SRR5206817: r14-resist-noninf (28 of 50)"
fastq-dump --stdout --gzip SRR5206817 > r14-resist-noninf.fastq.gz
echo "Downloading SRR5206818: r15-resist-infect (29 of 50)"
fastq-dump --stdout --gzip SRR5206818 > r15-resist-infect.fastq.gz
echo "Downloading SRR5206819: r15-resist-noninf (30 of 50)"
fastq-dump --stdout --gzip SRR5206819 > r15-resist-noninf.fastq.gz
echo "Downloading SRR5206820: r16-resist-infect (31 of 50)"
fastq-dump --stdout --gzip SRR5206820 > r16-resist-infect.fastq.gz
echo "Downloading SRR5206821: r16-resist-noninf (32 of 50)"
fastq-dump --stdout --gzip SRR5206821 > r16-resist-noninf.fastq.gz
echo "Downloading SRR5206822: r17-resist-infect (33 of 50)"
fastq-dump --stdout --gzip SRR5206822 > r17-resist-infect.fastq.gz
echo "Downloading SRR5206823: r17-resist-noninf (34 of 50)"
fastq-dump --stdout --gzip SRR5206823 > r17-resist-noninf.fastq.gz
echo "Downloading SRR5206824: r18-resist-infect (35 of 50)"
fastq-dump --stdout --gzip SRR5206824 > r18-resist-infect.fastq.gz
echo "Downloading SRR5206825: r18-resist-noninf (36 of 50)"
fastq-dump --stdout --gzip SRR5206825 > r18-resist-noninf.fastq.gz
echo "Downloading SRR5206826: r19-resist-infect (37 of 50)"
fastq-dump --stdout --gzip SRR5206826 > r19-resist-infect.fastq.gz
echo "Downloading SRR5206827: r19-resist-noninf (38 of 50)"
fastq-dump --stdout --gzip SRR5206827 > r19-resist-noninf.fastq.gz
echo "Downloading SRR5206828: s01-suscep-infect (39 of 50)"
fastq-dump --stdout --gzip SRR5206828 > s01-suscep-infect.fastq.gz
echo "Downloading SRR5206829: s01-suscep-noninf (40 of 50)"
fastq-dump --stdout --gzip SRR5206829 > s01-suscep-noninf.fastq.gz
echo "Downloading SRR5206830: s02-suscep-infect (41 of 50)"
fastq-dump --stdout --gzip SRR5206830 > s02-suscep-infect.fastq.gz
echo "Downloading SRR5206831: s02-suscep-noninf (42 of 50)"
fastq-dump --stdout --gzip SRR5206831 > s02-suscep-noninf.fastq.gz
echo "Downloading SRR5206832: s03-suscep-infect (43 of 50)"
fastq-dump --stdout --gzip SRR5206832 > s03-suscep-infect.fastq.gz
echo "Downloading SRR5206833: s03-suscep-noninf (44 of 50)"
fastq-dump --stdout --gzip SRR5206833 > s03-suscep-noninf.fastq.gz
echo "Downloading SRR5206834: s04-suscep-infect (45 of 50)"
fastq-dump --stdout --gzip SRR5206834 > s04-suscep-infect.fastq.gz
echo "Downloading SRR5206835: s04-suscep-noninf (46 of 50)"
fastq-dump --stdout --gzip SRR5206835 > s04-suscep-noninf.fastq.gz
echo "Downloading SRR5206836: s05-suscep-infect (47 of 50)"
fastq-dump --stdout --gzip SRR5206836 > s05-suscep-infect.fastq.gz
echo "Downloading SRR5206837: s05-suscep-noninf (48 of 50)"
fastq-dump --stdout --gzip SRR5206837 > s05-suscep-noninf.fastq.gz
echo "Downloading SRR5206838: s06-suscep-infect (49 of 50)"
fastq-dump --stdout --gzip SRR5206838 > s06-suscep-infect.fastq.gz
echo "Downloading SRR5206839: s06-suscep-noninf (50 of 50)"
fastq-dump --stdout --gzip SRR5206839 > s06-suscep-noninf.fastq.gz
