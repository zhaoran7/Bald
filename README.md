# Bald archaic-introgression haplotype analysis

本项目整理了男性秃发 GWAS/COJO 独立位点的古人类基因渗入候选片段分析流程：从 COJO lead SNP 出发，定义高 LD core risk haplotype，比较 1000G 与古人类单倍型，计算 ILS，并构建系统发育树。

## 数据准备

只需要准备 4 个 GWAS summary statistics 到 `data/gwas/`，并使用本仓库已包含的 `data/cojo/`；GWAS 文件名为 `bald.qt.gz`、`bald12.bt.gz`、`bald13.bt.gz`、`bald14.bt.gz`，至少包含 `SNP/CHR/POS/EA/NEA/BETA/P`。

```text
Bald/
├── data/
│   ├── cojo/      # 已包含 COJO 结果
│   └── gwas/      # 放入 4 个 GWAS 文件
├── scripts/       # s0-s8 分析脚本
└── README.md
```

脚本中保留了原 HPC 路径；如在其他环境运行，统一修改各脚本开头的 `proj`、`res`、`root` 等路径即可。

## 环境

需要 `bash`、`bcftools`、`tabix`、`plink`、`phyml`、`R`，以及 R 包 `data.table`、`ape`、`eoffice`。本流程默认使用 hg19 坐标。

## 流程

### s0. 下载参考数据

```bash
bash scripts/s0_download.sh
```

下载/检查 1000 Genomes phase 3、sample panel、Vindija、Altai、Chagyrskaya 和 Denisova VCF。

### s0_prepare. 准备 lead 与 risk allele

```bash
Rscript scripts/s0_prepare_inputs.R
```

从 `data/cojo/<trait>/jma.cojo` 和 `data/gwas/<trait>.gz` 生成：

```text
data/lead/<trait>.lead.tsv
res/risk.tsv
```

`res/risk.tsv` 中 risk allele 的定义为：`BETA >= 0` 取 effect allele，`BETA < 0` 取 non-effect allele。

### s1. 计算 lead SNP 的高 LD core region

```bash
bash scripts/s1_sumbit_ld.sh
bash scripts/s1_merge_ld.sh
```

每个 lead SNP 用 1000G VCF 计算 1 Mb 窗口内 LD：

```text
plink --r2 --ld-window-kb 1000 --ld-window-r2 0.98
```

输出：

```text
res/ld/<trait>/ld.tsv
res/ld/<trait>/block.tsv
```

`core region` 定义为与 lead SNP `r2 >= 0.98` 的 SNP 覆盖区间。

### s2. 保留多 SNP block

```bash
bash scripts/s2_pick.sh
```

筛选条件：

```text
n_ld_snp > 1
```

理由：单 SNP block 无法构建多位点 risk haplotype，也不能用于可靠的树分析。

输出：

```text
res/pick/pick.tsv
```

### s3. 提取 1000G 与古人类 core VCF

```bash
bash scripts/s3_submit_extract.sh
```

对每个 candidate block 提取 1000G、Vindija、Altai、Chagyr、Denisova 的 VCF：

```text
res/core_vcf/<trait>/<chr>.<lead_snp>.<lead_bp>/
```

### s4. VCF 转矩阵

```bash
bash scripts/s4_submit.sh
```

输出：

```text
res/mat/<trait>/<id>/kg.tsv
res/mat/<trait>/<id>/vindija.tsv
res/mat/<trait>/<id>/altai.tsv
res/mat/<trait>/<id>/chagyr.tsv
res/mat/<trait>/<id>/denisova.tsv
```

`kg.tsv` 保留 1000G phased genotype 和 ancestral allele；古人类矩阵保留 genotype。

### s5. 定义 risk haplotype 并筛选 archaic-matched loci

```bash
bash scripts/s5_submit.sh
```

核心逻辑在 `scripts/s5_make_hap.R`：

1. 对每个 lead SNP，在 `r2 >= 0.98` 的 core SNP 上提取 1000G phased haplotypes。
2. 找出携带 lead risk allele 的 haplotypes。
3. 对每个 core SNP，在 lead-risk haplotypes 中取多数等位基因，定义 `risk_core_allele`。
4. 将古人类 genotype 只在 homozygous 时转换为 allele；heterozygous/missing 记为 `NA`。
5. 比较古人类 allele 与 `risk_core_allele`。
6. 保留满足 `p_ils < 0.1` 且至少一个古人类样本通过匹配阈值的 region。

古人类匹配阈值：

```text
n_compared_risk >= 2
n_match_risk >= 2
prop_match_risk >= 0.5
```

ILS 使用固定参数估计：

```text
recombination rate = 0.53 cM/Mb
split time = 550,000 years
archaic age = 50,000 years
generation time = 29 years
```

输出：

```text
res/core_archaic_match.tsv
res/core_risk.tsv
res/hap_match.tsv
res/hap_site_count.tsv
res/region_summary.tsv
```

### s6. 生成 PHYLIP 输入

```bash
Rscript scripts/s6_make_phy.R
```

为每个候选 region 生成：

```text
res/phy/<trait>/<id>.full.phy
res/phy/<trait>/<id>.main.phy
```

`full` 包含所有保留 haplotypes；`main` 只保留出现次数 `n > 10` 的主要 haplotypes；同时加入古人类序列和 ancestral root。

### s7. PhyML 建树

```bash
bash scripts/s7_phyml_submit.sh
```

参数：

```text
phyml -m HKY85 -c 4 -a e -v e -b 100
```

输出：

```text
res/phy/<trait>/<id>.*.phy_phyml_tree.txt
```

### s8. 绘制树图

```bash
Rscript scripts/s8_plot_tree.R
```

输出：

```text
res/plot/s8_tree_full.pptx
res/plot/s8_tree_main.pptx
```

图中保留 PhyML 原始 branch length，并标注古人类样本、ancestral root 和高 bootstrap clade。

## 当前数据集的筛选结果

数量单位为 trait-region 记录；最后合并相同 genomic id 得到 13 个 region。

| 步骤 | 保留 | 筛掉 | 理由 |
|---|---:|---:|---|
| COJO independent lead signals | 1114 | - | 起始 lead SNP |
| Multi-SNP LD blocks | 555 | 559 | 去除 `n_ld_snp = 1` 的 singleton block |
| ILS filter | 117 | 438 | 保留 `p_ils < 0.1` |
| Archaic match filter | 60 | 57 | 去除 callable risk SNP 太少或匹配比例/数量不足的位点 |
| 可构建 haplotype/tree | 18 | 42 | 去除 informative sites 不足、古人类 callable 位点不足或无重复现代 haplotype 的记录 |
| 合并 genomic id | 13 | - | 多 trait 命中的同一 region 合并 |

最终 13 个 region 中：

```text
Neanderthal-like: 7
Denisovan-like:   6
```
