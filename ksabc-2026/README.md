# KSABC 2026 포스터 — 지루성 두피의 공존 네트워크와 *Lawsonella*

**Sebum-linked reorganization of scalp microbial co-occurrence networks and the rise of lipophilic *Lawsonella* in seborrheic dermatitis**

2026 한국응용생명화학회(KSABC) 포스터 발표의 그림 코드입니다.

> ## ⚠️ 이 포스터의 핵심 주장 하나는 철회했습니다
>
> 포스터의 헤드라인인 **"지루성 두피에서 공존 네트워크가 붕괴한다"(엣지 77 → 34, F3)는 이후 재검정을 통과하지 못했습니다.**
> 상대풍부도 자료에 맞는 방법(CLR 변환 + graphical lasso)으로 다시 계산하면 재현되지 않고, 방향이 오히려 뒤집혔습니다.
> 재검정 코드는 [`msk-2026/analysis/03_network_test.py`](../msk-2026/analysis/03_network_test.py) 에 있습니다.
>
> 같은 데이터의 후속 결론(피지·수분 두 독립 축)은 [`msk-2026/`](../msk-2026/) 를 보세요.
> 이 폴더는 **포스터를 그때 어떻게 만들었는지 기록으로 남기려고** 공개합니다.

## 폴더 구성

```
ksabc-2026/
├── scripts/                       포스터 그림 (번호 순서대로 실행)
│   ├── 00_prep_data.R               입력 표 → data/prep.rds (속 단위 표·메타데이터)
│   ├── 00_setup.R                   공통 색·테마·저장 함수 (모든 스크립트가 source)
│   ├── 01_F1_diversity_envfit.R     F1  알파 다양성(ACE·Shannon·Simpson) + Bray–Curtis PCoA·envfit
│   ├── 02_F2_anno_heatmap.R         F2  상위 25속 히트맵 (군 안에서 피지 순 정렬, 문(phylum)·시료 주석)
│   ├── 03_F3_diff_network.R         F3  NS / SS 공존 네트워크 비교  ← 철회 (위 참고)
│   ├── 03b_F3_themed.R              F3  같은 분석, 포스터용 디자인
│   ├── 04_F4_spearman.R             F4  상위 30속 × 임상 지표 Spearman 상관
│   ├── 05_F5_ridge.R                F5  피지 4분위별 주요 속 분포 (ridge plot)
│   ├── 06_F6_picrust2.R             F6  PICRUSt2 예측 대사경로 — 지질 관련 경로
│   └── install_*.R                  필요한 R 패키지 설치
└── v2_redesign/                   F1·F3·F6 를 다른 형식으로 다시 그린 시안
    ├── v2_F1a_raincloud.R           알파 다양성 → raincloud plot
    ├── v2_F3b_dumbbell.R            네트워크 지표(엣지 수·평균 degree·밀도·모듈성) NS vs SS → dumbbell plot
    └── v2_F6_lollipop.R             PICRUSt2 경로 → lollipop plot
```

모든 스크립트는 **`ksabc-2026/` 폴더에서** 실행합니다.

```bash
cd ksabc-2026
Rscript scripts/00_prep_data.R
Rscript scripts/01_F1_diversity_envfit.R
```

결과는 `output/`(시안은 `output/v2/`)에 저장됩니다.

## 해석할 때 주의할 점

- **F3 네트워크 결과는 철회했습니다** (위 상자 참고). Spearman |r| ≥ 0.6, p < 0.05 기준의 상관 네트워크는 상대풍부도 자료의 구성 효과(compositionality)에 취약합니다.
- **F4 상관 히트맵의 별표(\*)는 다중검정 보정 전 p 값 기준입니다.** 30속 × 임상 지표 여러 개를 동시에 검정했으므로, 별표 하나하나를 독립적인 발견으로 읽으면 안 됩니다.
- **F6 은 보정 p(BH)를 계산하지만 별표는 보정 전 p 값 기준입니다.** 또한 PICRUSt2 결과는 16S 로부터 **예측한** 기능이지 측정한 기능이 아닙니다.
- 표본은 27개(정상 두피 NS / 지루성 두피 SS)입니다.

## 데이터

**데이터는 이 저장소에 없습니다.** 사람 대상 시료라 공개하지 않습니다. 스크립트는 `ksabc-2026/data/` 에서 읽습니다.

| 입력 | 내용 |
|---|---|
| `data/feature_table_rare.tsv`, `data/feature_table_raw.tsv` | QIIME 2 에서 내보낸 ASV 표 (rarefied / raw) |
| `data/taxonomy/taxonomy.tsv` | ASV 분류 |
| `data/metadata_sebum.txt` | 군(NS/SS)·임상 지표 |
| `data/picrust2_out/` | PICRUSt2 결과 (`pathways_out/path_abun_unstrat.tsv.gz`) |

선택: 환경변수 `PICRUST2_METACYC_DESC` 에 PICRUSt2 의 `metacyc_pathways_info.txt.gz` 경로를 주면 F6 에 경로 설명이 붙습니다.

## 환경

R 패키지: ggplot2, vegan, reshape2, ggrepel, ape, igraph, patchwork, ggridges, ggnewscale, dplyr, tidyr, Hmisc — `scripts/install_*.R` 로 설치합니다(설치 위치 `~/R_libs`).
