# MSK 2026 포스터 — 두피 세균 군집의 두 축

**Two independent physiological axes structure the scalp bacterial community: hydration governs diversity, sebum governs the balance among dominant genera**

2026 International Meeting of the Microbiological Society of Korea (MSK 2026) 포스터 발표의 분석·그림 코드입니다.

## 요약

두피 시료 27개(정상 두피 NS / 지루성 두피 SS)의 16S rRNA 군집을 임상 지표와 함께 분석했습니다.

- **수분 축**: 수분이 군집 **다양성**을 결정한다 (나이 통제 후에도 유지).
- **피지 축**: 피지가 우점 3속 사이의 **균형**을 결정한다.
- 두 축은 서로 독립이다(피지 ~ 수분, 나이 통제 편상관 ≈ 0).
- 독립 코호트 2개(PRJNA1465550, PRJNA1199187)에서 사전 예측을 검증했다. 예측은 실패했고, 그 실패가 두 축 구조로 설명된다.

## 폴더 구성

```
msk-2026/
├── analysis/                  분석 (번호 순서대로 실행)
│   ├── 01_two_axis.py           두 축 분리 — 확정 수치
│   ├── 02_varpart.R             분산분할 — 피지·수분 각각의 고유 설명력
│   ├── 03_network_test.py       KSABC 포스터의 '공존 네트워크 붕괴' 주장 재검정 → 재현 안 됨
│   ├── 04_ml_axes.py            군집 조성으로 피지·수분 예측
│   ├── 05_external_validation.py  독립 코호트 외부 검증
│   ├── 06_server_ml_external.py   04+05 통합 (순열 999회 — 서버에서 실행)
│   ├── 07_reinterpret.py        외부 검증 실패의 해석 (사후 해석이므로 검정을 붙임)
│   ├── 08_ext2_PRJNA1199187/    두 번째 외부 코호트 (16S·ITS 직교성)
│   ├── 10_ext1_recheck/         첫 외부 코호트 재검 (baseline/treatment 분리, 짝 검정)
│   ├── 12_sample_structure/     유사반복(pseudoreplication) 검정 — 27개 시료가 27명인가
│   ├── 14_all_indices/          임상 지표 7개 전부 검정 — 정말 두 축뿐인가
│   └── 17_small_n/              작은 표본 보강 — 효과크기 + 95% 신뢰구간, 포함률 시뮬레이션
└── figures/                   포스터 그림 F1–F6
    ├── fig_style.py             공통 스타일·데이터 로더
    ├── F1_two_axes.py … F6_evidence.py
    └── F_legends.py             패널 범례
```

각 스크립트 맨 위 docstring 에 입력·출력·검정 원칙이 적혀 있습니다.

## 데이터

**데이터는 이 저장소에 없습니다.** 사람 대상 시료라 원자료와 임상 지표는 공개하지 않습니다.
스크립트는 `msk-2026/data/` 에서 입력을 읽고 `msk-2026/output/` 에 결과를 씁니다.

| 입력 | 내용 | 공개 여부 |
|---|---|---|
| `data/feature-table.biom`, `data/taxonomy.tsv`, `data/metadata_sebum.txt`, `data/own_asv_table.tsv` | 자체 코호트 ASV 표·분류·임상 지표 | 비공개 |
| `data/external/own_taxonomy_silva138.tsv` | 자체 코호트를 외부 코호트와 같은 분류기(SILVA 138)로 재분류한 것 | 비공개 |
| `data/external/feature-table.tsv`, `data/external/taxonomy.tsv` | 외부 코호트 1 (PRJNA1465550) 의 ASV 표·분류 | NCBI SRA 공개 |
| `data/PRJNA1465550/` | 외부 코호트 1 — QIIME2 `export_table/`·`export_taxonomy/` (10_ext1_recheck 용) | NCBI SRA 공개 |
| `data/PRJNA1199187/` | 외부 코호트 2 — 논문 보충표 (Table S1) | 논문 공개 |
| `analysis/10_ext1_recheck/runs_*.tsv` | 외부 코호트 1 의 run ↔ 대상자·병변 대응표 | SRA 메타데이터에서 정리 (포함) |

## 환경

```bash
pip install numpy pandas scipy statsmodels scikit-learn biom-format matplotlib networkx joblib
# R: install.packages("vegan")
```

난수 시드는 전부 `20260825` 로 고정했습니다.
