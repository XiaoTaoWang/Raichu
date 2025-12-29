Raichu 
======
Accurate detection of enhancer-promoter loops from genome-wide chromatin interaction
data is critical for understanding gene regulation. Standard normalization methods,
such as matrix balancing approaches, are widely used to correct biases in chromatin
contact data prior to chromatin loop detection. However, while these methods preserve
structural loop signals, they often attenuate enhancer-promoter interaction signals,
making these regulatory loops more difficult to detect. To address this limitation,
we develop Raichu, a novel normalization method for chromatin contact data. Raichu
identifies nearly twice as many loops as conventional normalization approaches,
recovering almost all previously detected loops while uncovering thousands of
additional enhancer-promoter interactions that are otherwise missed. With its
improved sensitivity for regulatory loops, Raichu detects more biologically meaningful
differential interactions, including those between conditions within the same cell
type. Moreover, Raichu performs robustly across a wide range of sequencing depths,
resolutions, species, and experimental platforms, making it a versatile tool for
revealing new insights into three-dimensional genome organization and transcriptional
regulation.

Citation
========
Wang, X., Shi, D., Xue, F., Liu, Y., Yang, H., Jiang, L. Boosting the detection of
enhancer-promoter loops via normalization methods for chromatin interaction data.
Nature Communications. 2026.

Installation
============
Raichu and all the dependencies can be installed through either `mamba <https://github.com/mamba-org/mamba>`_
or `pip <https://pypi.org/project/pip/>`_::

    $ conda config --append channels defaults
    $ conda config --append channels bioconda
    $ conda config --append channels conda-forge
    $ mamba create -n 3Dnorm numba joblib "cooler==0.9.3" "scipy>=1.10"
    $ mamba activate 3Dnorm
    $ pip install RaichuNorm

Raichu is a command-line tool, and after successful installation, help information
can be accessed by running ``raichu -h`` in a terminal.

Usage
=====
Raichu is built on the `cooler <https://github.com/open2c/cooler>`_ Python package
for reading and processing contact matrices. To demonstrate how to normalize a
contact matrix in .cool format, let's download the file "GM12878.Hi-C.10kb.cool"
from this `link <https://www.jianguoyun.com/p/DUoSz7gQh9qdDBi5lLwFIAA>`_. This
file contains contact matrices at 10kb resolution, generated from an in situ Hi-C
dataset in the GM12878 cell line.

Now all that is needed is to execute the commands below in a terminal::

    $ raichu --cool-uri GM12878.Hi-C.10kb.cool --window-size 200 -p 8 -n obj_weight -f

Here:

1. The ``--cool-uri`` parameter specifies the URI of contact matrices at
a specific resolution. For a single-resolution cooler file (typically suffixed
with .cool), the value should be the file path. For a multi-resolution cooler
file (typically suffixed with .mcool), the value should include the file path
followed by ``::`` and the internal group path to the root of a data collection.
For example: ``test.mcool::resolutions/10000`` or ``test.mcool::resolutions/5000``.

2. The ``--window-size`` parameter specifies the size of the sliding window. In most
cases, the default value of 200 is sufficient. Increasing the window size may
improve the accuracy of bias vector calculations but will also increase the runtime.

3. The ``-p`` or ``--nproc`` parameter specifies the number of processes to allocate for
the calculation. Raichu uses this parameter to perform calculations for chromosomes
in parallel. However, setting this parameter to a value greater than the number of
chromosomes will not result in additional speed improvements.

4. The ``-n`` or ``--name`` parameter specifies the name of the column where the
calculated bias vectors will be written.

5. If the ``-f`` or ``--force`` parameter is specified, the target column in the
bin table will be overwritten if it already exists.

Downstream Analysis with Raichu-Normalized Matrices
===================================================
Raichu stores the calculated bias vectors in the same format as
``cooler balance`` (an implementation of the ICE algorithm), ensuring
seamless compability with downstream tools for analyzing compartments,
TADs, and loops.

For instance, to compute chromatin compartment values based on Raichu-normalized
signals, we can use the `cooltools eigs-cis  <https://github.com/open2c/cooltools>`_
command and specify the ``--clr-weight-name`` parameter as "obj_weight" (matching
the ``-n`` parameter setting we used when running Raichu). The full command would
look like this::

    $ cooltools eigs-cis --phasing-track hg38-gene-density-100K.bedGraph --clr-weight-name obj_weight -o GM_raichu GM12878-MboI-allReps-hg38.mcool::resolutions/100000

Similarly, we can use the following command to compute insulation scores with
Raichu-normalized signals::

    $ cooltools insulation --ignore-diags 1 -p 8 -o GM_raichu.IS.25kb.tsv --clr-weight-name obj_weight GM12878-MboI-allReps-hg38.mcool::resolutions/25000 1000000

For loop detection, we have tested the `pyHICCUPS <https://github.com/XiaoTaoWang/HiCPeaks>`_,
`Mustache <https://github.com/ay-lab/mustache>`_, and `Peakachu <https://github.com/tariks/peakachu>`_
software. See the next section for the specific commands used in our analyses.

Loop-calling commands used in the manuscript
============================================
This section documents the exact commands used to identify chromatin loops in the datasets
analyzed in the manuscript, enabling full reproducibility of the benchmarking and comparative
analyses. Note that for each dataset, the same parameters were applied when using ICE-normalized
signals for loop detection, except that the ``--clr-weight-name`` parameter was set to "weight".

GM12878 Hi-C data (Figure 2)
----------------------------
To identify loops from the GM12878 Hi-C data, we used the following commands with pyHICCUPS v0.3.9::

    $ pyHICCUPS -p GM12878-MboI-allReps-hg38.mcool::resolutions/5000 -O GM12878_pyHICCUPS.5kb.bedpe \
                --pw 1 2 4 --ww 3 5 7 --only-anchors --nproc 8 --clr-weight-name obj_weight --maxapart 4000000
    $ pyHICCUPS -p GM12878-MboI-allReps-hg38.mcool::resolutions/10000 -O GM12878_pyHICCUPS.10kb.bedpe \
                --pw 1 2 4 --ww 3 5 7 --only-anchors --nproc 8 --clr-weight-name obj_weight --maxapart 4000000
    $ combine-resolutions -O GM12878_pyHICCUPS.bedpe -p GM12878_pyHICCUPS.5kb.bedpe GM12878_pyHICCUPS.10kb.bedpe -R 5000 10000 -G 10000 -M 100000 --max-res 10000

GM12878 Hi-C (Mustache)
-----------------------
Loops were also identified from the same GM12878 Hi-C dataset using Mustache v1.3.2
(Supplementary Figure 14)::

    $ mustache -f GM12878-MboI-allReps-hg38.mcool -r 5000 -pt 0.05 -norm obj_weight -p 8 -o GM12878_mustache.5kb.tsv
    $ mustache -f GM12878-MboI-allReps-hg38.mcool -r 10000 -pt 0.05 -norm obj_weight -p 8 -o GM12878_mustache.10kb.tsv
    $ combine-resolutions -O GM12878_mustache.bedpe -p GM12878_mustache.5kb.tsv GM12878_mustache.10kb.tsv -R 5000 10000 -G 10000 -M 100000 --max-res 10000

mESC Hi-C (Figure 3)
--------------------
To identify loops from the mESC Hi-C data, we used pyHICCUPS v0.3.9 with the following commands::

    $ pyHICCUPS -p HIC-mESC-DpnII-allReps.10kb.cool -O HIC-mESC-hicpeaks.10kb.raichu.min1000.txt \
                --pw 2 4 --ww 5 7 --maxww 7 --only-anchors --min-local-reads 1000 --min-marginal-peaks 3 \
                --nproc 8 --clr-weight-name obj_weight --maxapart 3000000 --logFile HIC-mESC-hicpeaks.raichu.log
    $ pyHICCUPS -p HIC-mESC-DpnII-allReps.5kb.cool -O HIC-mESC-hicpeaks.5kb.raichu.min1000.txt \
                --pw 2 4 --ww 5 7 --maxww 7 --only-anchors --min-local-reads 1000 --min-marginal-peaks 3 \
                --nproc 8 --clr-weight-name obj_weight --maxapart 1000000 --logFile HIC-mESC-hicpeaks.raichu.log
    $ combine-resolutions -O HIC-mESC-hicpeaks.raichu.min1000.bedpe -p HIC-mESC-hicpeaks.5kb.raichu.min1000.txt HIC-mESC-hicpeaks.10kb.raichu.min1000.txt \
                          -R 5000 10000 -G 10000 -M 100000 --max-res 10000

mNPC and hNPC Hi-C (Figure 6)
-----------------------------
To identify loops from the mNPC Hi-C data, we used pyHICCUPS v0.3.9 with the following commands::

    $ pyHICCUPS -p HIC-mNPC-DpnII-allReps.5kb.cool -O HIC-mNPC-hicpeaks.5kb.raichu.min300.txt \
                --pw 1 2 4 --ww 3 5 7 --only-anchors --min-local-reads 300 --nproc 8 --clr-weight-name obj_weight \
                --maxapart 4000000 --logFile HIC-mNPC-hicpeaks.raichu.log
    $ pyHICCUPS -p HIC-mNPC-DpnII-allReps.10kb.cool -O HIC-mNPC-hicpeaks.10kb.raichu.min300.txt \
                --pw 1 2 4 --ww 3 5 7 --only-anchors --min-local-reads 300 --nproc 8 --clr-weight-name obj_weight \
                --maxapart 4000000 --logFile HIC-mNPC-hicpeaks.raichu.log
    $ combine-resolutions -O mNPC-raichu.min300.combined.bedpe -p HIC-mNPC-hicpeaks.5kb.raichu.min300.txt HIC-mNPC-hicpeaks.10kb.raichu.min300.txt \
                          -R 5000 10000 -G 10000 -M 100000 --max-res 10000

The same parameters were applied to the hNPC Hi-C data.

H1ESC Micro-C (Figure 7)
------------------------
Loops from the H1ESC Micro-C data were identified using pyHICCUPS v0.3.9 as well::

    $ pyHICCUPS -p MicroC-H1ESC-MNase-allReps.10kb.cool -O H1ESC_MicroC_hicpeaks.10kb.raichu.min100.txt \
                --pw 2 3 4 --ww 5 6 7 --maxww 7 --only-anchors --min-local-reads 100 --nproc 8 \
                --clr-weight-name obj_weight --maxapart 4000000 --logFile H1ESC_MicroC_hicpeaks.raichu.log
    $ pyHICCUPS -p MicroC-H1ESC-MNase-allReps.5kb.cool -O H1ESC_MicroC_hicpeaks.5kb.raichu.min100.txt \
                --pw 2 3 4 --ww 5 6 7 --maxww 7 --only-anchors --min-local-reads 100 --nproc 8 \
                --clr-weight-name obj_weight --maxapart 2000000 --logFile H1ESC_MicroC_hicpeaks.raichu.log
    $ pyHICCUPS -p MicroC-H1ESC-MNase-allReps.2kb.cool -O H1ESC_MicroC_hicpeaks.2kb.raichu.min100.txt \
                --pw 2 3 4 --ww 5 6 7 --maxww 7 --only-anchors --min-local-reads 100 --nproc 8 \
                --clr-weight-name obj_weight --maxapart 1000000 --logFile H1ESC_MicroC_hicpeaks.raichu.log
    $ combine-resolutions -O H1ESC_MicroC_hicpeaks.raichu.min100.bedpe -p H1ESC_MicroC_hicpeaks.2kb.raichu.min100.txt \
                          H1ESC_MicroC_hicpeaks.5kb.raichu.min100.txt H1ESC_MicroC_hicpeaks.10kb.raichu.min100.txt \
                          -R 2000 5000 10000 -G 10000 -M 100000 --max-res 10000 --minapart 8

mESC region-capture Micro-C (Figure 7)
--------------------------------------
First, we constructed a BED file defining the captured regions::

    $ cat capture-regions.bed

    chr5    31257344        32382344
    chr8    84846629        85856629
    chr18   58032072        59034072
    chr3    33804149        35704149
    chr6    122451959       122876959

Bias vectors were then calculated only within the captured regions using Raichu::

    $ raichu --cool-uri GSM6281849_RCMC_BR1_merged_allCap_WT_mm39.merged.50.mcool::resolutions/1000 \
             --window-size 800 --regions capture-regions.bed --ignore-diags 0 --lower-bound 0.001 \
             --upper-bound 1000 -p 4 -n obj_weight -f --logFile RCMC_raichu.log
    $ raichu --cool-uri GSM6281849_RCMC_BR1_merged_allCap_WT_mm39.merged.50.mcool::resolutions/500 \
             --window-size 800 --regions capture-regions.bed --ignore-diags 0 --lower-bound 0.001 \
             --upper-bound 1000 -p 4 -n obj_weight -f --logFile RCMC_raichu.log
    $ raichu --cool-uri GSM6281849_RCMC_BR1_merged_allCap_WT_mm39.merged.50.mcool::resolutions/250 \
             --window-size 800 --regions capture-regions.bed --ignore-diags 0 --lower-bound 0.001 \
             --upper-bound 1000 -p 4 -n obj_weight -f --logFile RCMC_raichu.log

Finally, loops within the captured regions were identified using pyHICCUPS v0.3.9::

    $ pyHICCUPS -p GSM6281849_RCMC_BR1_merged_allCap_WT_mm39.merged.50.mcool::resolutions/1000 \
                -O mESC_RCMC_hicpeaks.1kb.raichu.min100.txt --pw 2 3 4 --ww 5 6 7 --maxww 7 \
                -C 5 8 18 3 6 --only-anchors --min-local-reads 100 --nproc 1 --clr-weight-name obj_weight \
                --maxapart 400000 --logFile mESC_RCMC_hicpeaks.log
    $ pyHICCUPS -p GSM6281849_RCMC_BR1_merged_allCap_WT_mm39.merged.50.mcool::resolutions/500 \
                -O mESC_RCMC_hicpeaks.500bp.raichu.min100.txt --pw 2 3 4 --ww 5 6 7 --maxww 7 \
                -C 5 8 18 3 6 --only-anchors --min-local-reads 100 --nproc 1 --clr-weight-name obj_weight \
                --maxapart 200000 --logFile mESC_RCMC_hicpeaks.log
    $ pyHICCUPS -p GSM6281849_RCMC_BR1_merged_allCap_WT_mm39.merged.50.mcool::resolutions/250 \
                -O mESC_RCMC_hicpeaks.250bp.raichu.min100.txt --pw 2 3 4 --ww 5 6 7 --maxww 7 \
                -C 5 8 18 3 6 --only-anchors --min-local-reads 100 --nproc 1 --clr-weight-name obj_weight \
                --maxapart 100000 --logFile mESC_RCMC_hicpeaks.log
    $ combine-resolutions -O mESC_RCMC_hicpeaks.raichu.min100.bedpe -p mESC_RCMC_hicpeaks.250bp.raichu.min100.txt \
                          mESC_RCMC_hicpeaks.500bp.raichu.min100.txt mESC_RCMC_hicpeaks.1kb.raichu.min100.txt \
                          -R 250 500 1000 -G 500 -M 5000 --max-res 500 --minapart 8

Performance
===========
In GM12878 cells, ICE detected 15,446 loops, while Raichu identified 28,986 loops.
(For this analysis, pyHICCUPS was applied; however, as shown in the manuscript,
various loop-calling methods achieve a similar level of improvement when using
Raichu-normalized signals.) Notably, 90.6% of loops detected by ICE (13,997 out
of 15,446) were also identified by Raichu, whereas 51.7% of loops detected by
Raichu (14,989 out of 28,986) were missed by ICE.

We classified the loops into three categories: ICE-specific loops, Raichu-specific loops,
and common loops (detected by both ICE and Raichu). While ICE-specific, Raichu-specific,
and common loops showed comparable enrichment for CTCF and RAD21, Raichu-specific
loops exhibited substantially greater enrichment for a broader range of transcription
factors (TFs) and histone modifications closely associated with transcriptional regulation.

.. image:: ./images/performance.png
        :align: center
