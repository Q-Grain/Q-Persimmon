# gRNA 候选设计 v0.1

- 日期：2026-10-28
- 底盘：君迁子 *Diospyros lotus*（`GCF_014633365.1`，2n=30）
- 脚本：`scripts/design_grna_candidates.ps1`、`scripts/locate_guide_hits.ps1`
- 产物：`data/targets/grna_candidates.tsv`（52 条）、`grna_offtarget.tsv`（带脱靶列）、`grna_shortlist.tsv`（15 条短名单）

## 1. 靶点怎么选的

从 `target-genes.md` 的 36 个核心里，只用满足以下条件的（**只选"敲掉有意义"的降解侧**）：

| 条件 | 原因 |
|---|---|
| 是淀粉**降解**酶 | 增强侧要靠过表达/启动子编辑，敲除没意义 |
| 注释为**活跃**酶 | `inactive beta-amylase 4 / 9`、`inactive AGPase 小亚基` 被排除——敲掉无活性的酶等于白干 |
| 定位在**质体/淀粉体** | 果实里的淀粉存在淀粉体中，酶得进得去才拆得动 |
| 模型**完整** | 排除 `LOC127789021`（342 bp）这类残缺模型 |
| 注释明确 | 排除 `probable`、`like`、`uncharacterized` 等不确定项 |

最终 5 个：

| 基因 | 角色 | 敲哪 | CDS 长度 |
|---|---|---|---|
| `LOC127799073` | BAM1（质体 β-淀粉酶 1） | 前 2 个外显子 | 1,734 |
| `LOC127786821` | BAM3（质体 β-淀粉酶 3） | 前 2 个外显子 | 1,638 |
| `LOC127807677` | ISA3（异淀粉酶 3，脱支） | 前 2 个外显子 | 2,334 |
| `LOC127803838` | PHS-L1（质体淀粉磷酸化酶） | 前 2 个外显子 | 2,994 |
| `LOC127806752` | PUL1（普鲁兰酶 1，脱支） | 前 2 个外显子 | 2,907 |

**为什么打前 2 个外显子**：敲除的目的是让蛋白失活，打在前段外显子最容易造成移码、把后面的结构域全废掉。

## 2. 短名单（15 条）

过滤条件：GC 40~70%、无 `TTTT`（Pol III 终止）、无 `GGGG`（G4）、无 BsaI/BsmBI 位点；同一基因内起始位置间隔 ≥10 nt。

| 基因 | 外显子 | 染色体 | 起点 | 链 | gRNA (20 nt) | PAM | GC% | 完全匹配 | 错配≤3 |
|---|---|---|---|---|---|---|---|---|---|
| ISA3 | 2 | chr8 | 2,134,629 | + | `TCACTCGGCCCTTTCTTGCT` | CGG | 55 | 1 | 1 |
| ISA3 | 1 | chr8 | 2,135,247 | + | `TTCATCGCCGAATCGCGAAG` | AGG | 55 | 1 | 1 |
| ISA3 | 2 | chr8 | 2,134,614 | + | `CTTTGATCTGTGTTGTCACT` | CGG | 40 | 1 | 1 |
| PUL1 | 1 | chr7 | 22,049,120 | + | `AGATGAGGCGTGGCAATGAA` | CGG | 50 | 1 | 1 |
| PUL1 | 1 | chr7 | 22,049,166 | + | `GCTTGTCTTTGGGAGGATGT` | CGG | 50 | 1 | 1 |
| PUL1 | 1 | chr7 | 22,049,217 | − | `CGATTCGTTACTTGACTCGC` | TGG | 50 | 1 | 1 |
| PHS-L1 | 1 | chr6 | 3,938,338 | + | `CGTTTCCTGTGACTTCGACT` | CGG | 50 | 1 | 1 |
| PHS-L1 | 1 | chr6 | 3,938,429 | − | `CTCTCGGGTCAAGAATAGCT` | TGG | 50 | 1 | 1 |
| PHS-L1 | 1 | chr6 | 3,938,439 | + | `TTGACCCGAGAGTCGAGTTT` | GGG | 50 | 1 | 1 |
| BAM3 | 1 | chr12 | 9,882,035 | + | `CCACTCGACACCATGTCATT` | TGG | 50 | 1 | 1 |
| BAM3 | 1 | chr12 | 9,882,036 | − | `CCAAATGACATGGTGTCGAG` | TGG | 50 | 1 | 1 |
| BAM3 | 1 | chr12 | 9,882,075 | − | `AAGCTGGCATTCATTGCTCG` | TGG | 50 | 1 | 1 |
| ⚠️ BAM1 | 1 | chr4 | 1,637,314 | + | `GTTTCCGCCACATTGATGGA` | AGG | 50 | **2** | 2 |
| ⚠️ BAM1 | 1 | chr4 | 1,637,826 | − | `GTGGATTCCGGTAATAACGC` | GGG | 50 | **2** | 2 |
| ⚠️ BAM1 | 1 | chr4 | 1,637,838 | − | `ATGCCGTTAGCGGTGGATTC` | CGG | 55 | **2** | 2 |

## 3. 重要发现：BAM1 那一段不能用

三条 BAM1 候选的"完全匹配"都是 2，不是 1。查下去发现：

| 检查 | 结果 |
|---|---|
| `GTTTCCGCCACATTGATGGA` 在基因组里的位置 | chr4:**1,637,314** 和 chr4:**1,661,360** |
| 两点距离 | 约 **24 kb** |
| 第二个点落在哪 | `LOC127799072` 的范围内 |
| 该区域结构 | `LOC127799071 / 072 / 073(BAM1) / 075` 依次排列，是一个**串联基因簇** |

**后果**：这条 gRNA 会同时在相距 24 kb 的两处切断 DNA，产生的不是一个小 indel，而是**跨 24 kb 的缺失或倒位**，细胞很可能承受不住，或者得到无法解释的复杂重排。

**结论**：BAM1 这三条**不推荐使用**。要做 BAM1，得回到基因组里重新找位于基因簇外侧、唯一性好的位点——这需要单独做一轮，不在这次范围内。

## 4. 方法说明与局限（必须知道）

| 项目 | 这次怎么做的 | 局限 |
|---|---|---|
| 候选生成 | 在前 2 个外显子内滑动找 20 nt + NGG，两条链都找 | —— |
| 脱靶搜索 | 8 nt 种子（PAM 侧）+ 全长 20 nt 比对 + PAM 必须匹配，允许 ≤3 错配 | 不支持 **bulge**（插入/缺失型错配），可能低估脱靶 |
| 脱靶评分 | 只统计位点**数量** | **没有 MIT / CFD 打分**，无法区分"错在种子区"和"错在远端" |
| 大小写 | 基因组有软屏蔽区（小写）。第一版漏筛了，已修 | —— |
| 靶点唯一性 | 用"完全匹配位点数"判断 | 只能发现完全相同的重复，发现不了高相似的旁系同源 |

一句话：**这批结果够用来筛，不够用来定稿。**

## 5. 下单前的复核清单

在合成之前，至少要过这四关：

1. **用 CHOPCHOP / CRISPR-P 复核这 12 条**（去掉 BAM1）的序列与位置，看它们的评分
2. **BLAST 一次**：确认没有高相似的旁系同源（我这边只查了完全匹配）
3. **用真实材料的 DNA 验证靶点**：从盘山柿/君迁子提 DNA，PCR 扩增靶区并测序，确认本地序列与参考基因组一致——**参考基因组和你的树可能不一样**
4. **确认 PAM 与载体匹配**：选 Cas9 变体（SpCas9 用 NGG）与 sgRNA 表达载体

## 6. 复现

```powershell
# 生成候选 + 全基因组脱靶 + 短名单（约 2~3 分钟）
powershell -NoProfile -ExecutionPolicy Bypass -File scripts\design_grna_candidates.ps1 -Offtarget

# 查某条 gRNA 在基因组里的所有完全匹配位置
powershell -NoProfile -ExecutionPolicy Bypass -File scripts\locate_guide_hits.ps1 -Guide GTTTCCGCCACATTGATGGA
```
