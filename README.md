# 🌊 FluidFullSpeedSimulation

💧 水・🍶 エタノール・🧴 メタノールなどの作動流体で、**流れ・遅れ・検証**を整数演算で回すシミュレーター。

- 🏛️ 手本は MONIAC（フィリップス・マシン, 1949）の「水が流れ落ちる様子」。経済の解釈は扱わず、部品（水槽・管・堰・浮き弁・遅延管・連結）を物理のまま再現する
- 🎋 部品の組み合わせで、流しそうめんの台のような別の題材にも広げられる汎用の回路にしている
- 🧮 量はすべて Q 整数（固定小数点）。移した量をそのまま足し引きするので、総量は 1 ulp も変わらない
- 🟦 見せ方は絵文字のブロックを積む断面図（行 = 高さ、列 = 区間・水槽）。溜まる・流れる・遅れるがそのまま見える

## 🚀 使い方

```bash
uv sync --locked                                       # 🔒 uv.lock どおり（実行は標準ライブラリだけ）
uv run python fluidfullspeedsimulation.py verify        # ✅ 全実験 → 終了コード 0 / 1
uv run python fluidfullspeedsimulation.py list          # 🗂️ 実験の一覧
uv run python fluidfullspeedsimulation.py fluids        # 🧪 流体の物性
uv run python fluidfullspeedsimulation.py run flume -v --jsonl plot.jsonl   # 🎋 1 本だけ + 記録
uv run python fluidfullspeedsimulation.py run 🥢                             # 実験はアイコンでも選べる
uv run python fluidfullspeedsimulation.py draw 🏛️ 🟩                         # 場面と流体（ブロック）を絵文字のまま
uv run python fluidfullspeedsimulation.py draw 🥢 🟥                          # 流体でない絵文字は、その色で塗るだけ
uv run python fluidfullspeedsimulation.py draw 🧊 --temp-c 5                  # 🧊 水温を変えて細い流れ（場面: 🏛️ 🥢 🧊 🪵 🎋 🫙 🌊 🕳️）
uv run python -m pytest -q                              # 🧪 テスト（python -S -m unittest discover -s tests でも可）
```

終了コード: `0` ✅ 全検証が通った / `1` ❌ どれかが外れた / `2` 引数の誤り。

### 🐍 / 🟩 名前で呼ぶ（`__main__` を通らない）

```python
from fluidfullspeedsimulation import run_experiment
run_experiment("E05", gain_high=0.03)          # 既定値を上書きして回す

from network import simulate                    # 🧩 回路を JSON の仕様で組む
simulate({"dt": 0.02, "tanks": [{"name": "sea", "area": 1, "volume": 1}],
          "channels": [{"prefix": "flume", "dst": "sea", "length": 2, "cells": 20, "slope": 0.05,
                        "section": {"shape": "semicircle", "size": 0.04}}],
          "links": [{"law": "pump", "src": "sea", "dst": "flume:0", "rate": 5e-5}],
          "tracers": [{"name": "somen", "node": "flume:0", "release": 15, "stop": "sea"}],
          "t_end": 25})
```

```bash
# vendor/py_launch.cjs（親リポジトリの py_launch.cjs を vendor.lock.json で固定。Node 26+）
node vendor/py_launch.cjs --call run_experiment -- fluidfullspeedsimulation.py E08
node vendor/py_launch.cjs --call simulate -- network.py '{"dt":1,"tanks":[...],"links":[...],"t_end":10}'
```

### 🔒 vendor

| ファイル | 元 | 使いどころ |
|---|---|---|
| `vendor/cli_args.py` | myon-bioinformatics/cli_args | CLI の parser（`make_parser`、省略形を受け付けない、終了コード 0/1/2） |
| `vendor/py_launch.cjs` | myon-bioinformatics/myon-bioinformatics | Node から関数を名前で呼ぶ |

commit・blob・sha256 は `vendor.lock.json` に固定。更新は親の `vendor_sync.py`（`check` / `update` / `enroll`）で。

戻り値は JSON 越しでも変わらない形（2^53 を超える整数は文字列、∞ は `"inf"`）。

## 🧩 部品

| | 部品 | 中身 |
|---|---|---|
| 🛢️ | `tank` | 水槽。水位 h = V / A。流入 − 流出を積み上げる積分器 |
| 🎋 | `channel` | 傾いた樋を区間セルに分けたもの。断面は長方形・半円（竹）。Manning（乱流）/ 層流の薄い流れ |
| 🥢 | `flume` | 竹を継いだ流しそうめんの台。継ぎ目は段差 `drop`（自由落下）か縁 `lip`（たまってから越える）、`bend_radius` の曲がりは外側が持ち上がって縁からあふれる |
| 〰️ | `pipe` | 細管（Hagen–Poiseuille）。速さは動粘度 ν で決まる |
| 🕳️ | `orifice` | 小孔（Torricelli）。流体によらない |
| 🪜 | `weir` / `cam` | 長方形・三角の堰、水位 → 流量の表（特殊形状のスリット） |
| ⛽ | `pump` | 一定流量（`until` まで） |
| 🎈 | `float_valve` | 浮き弁。測った水位が遅れの段を通って弁を開く |
| 🚰 | `valve` | 弁つきの出口（MONIAC の弁）。開き具合は一定か、別の水槽の浮き → カムで決まる |
| ⏳ | `delay` | どの管にも付く輸送の遅れ。管の中の水も総量に数える |
| 🍜 | `tracer` | 流れに乗って運ばれる物（そうめん・葉っぱ）。先入れ先出しで節点を渡る。太さ `size` があると、縁の越流水深が浅いあいだ引っかかる |
| 🔌 | `register_law` | 管の法則を足す口。spec の `law` からもそのまま使える |

## 🧪 実験（ID + 識別子）

| ID | 識別子 | 確かめること |
|---|---|---|
| E01 🧮 | `tank-integrator` | 一定流入 → 体積が整数で厳密に n·(Q·dt)、排水が e^(−t/τ) どおり、2^−s の遅れ = シフト |
| E02 💧 | `fluid-drain` | 細管は ν に比例して遅い（🍶 は 💧 の約 1.51 倍）、小孔は流体によらずビット単位で同じ |
| E03 ⏳ | `transport-delay` | 管を抜ける遅れ L/v と、管の中の水 |
| E04 🪜 | `weir` | 堰・カムの定常水位が閉じた式どおり、粘度によらない |
| E05 🎈 | `float-valve` | 浮き弁の遅れ 1 段は強い弁ほど設定値へ近く静か、3 段では限界ゲインを超えると揺れが残る |
| E06 🌐 | `coupled-loop` | 2 台を遅延管で環につなぐ: 定常の配分、届く時刻、跳ね返りが戻る時刻 |
| E07 🏎️ | `fluid-speed` | 同じ回路を流体だけ変える: 流量は同じ、水位と整定時間は ν の比 |
| E08 🎋 | `flume` | 樋の等流水深、トレーサー（そうめん役）が抜ける時間 L/v、層流は ν に効き Manning は効かない |
| E09 🟦 | `pixel-water` | 断面図: 段数から水深を ±½ 段で読み戻せる、先端の届く時刻の間隔 ≈ Δx/v、じょうごに溜まって段が積み上がる、🍜 は管の遅れの間 ⏳。ドット絵の水: ドット数が毎ティック不変、水面が平ら、静まる順 🧴 < 💧 < 🍶、seed で再現、PNG/GIF の読み戻し |
| E10 🥢 | `somen-flume` | 竹 3 本の台: 継ぎ目のたまり・段差の落下・流しすぎで床へあふれる・細い流れでひやむぎだけ引っかかりそうめんに抜かれる・水を足すと流れ出す |
| E11 🏛️ | `moniac-rig` | MONIAC 型の水回りを据え置きの試験台に: 3 つの弁の分流比、浮き → カムの弁で決まる定常水位、閉じている間は入った量がそのまま貯まる、別の槽の浮きで弁が開く、戻りの管の遅れ。最後の水位は `tests/test_moniac_rig.py` の GOLDEN と比べる |
| E15 🧊 | `cold-water` | 水温表（0–40 °C）: 出典間の幅、10 °C の細い流れは等流水深の比どおり深く（×1.08）遅い（×1.12）、縁の越流水深と引っかかりは変わらない、Manning の竹はビット単位で同じ |
| E16 🪵 | `bamboo-wear` | 竹の傷み: 割れ目の漏れ = Cd·a·√(2gh)、床に落ちたのは漏れた分だけ（ulp まで）、節の削り残しの手前のたまり = 節 + 越流水深、細い流れではうどん（2 mm）だけ節で引っかかる |

どの実験も ⚖️ 総量の保存（ずれ 0 ulp）と、層流の前提（Re）・Courant 数などの前提条件を一緒に確かめる。
閉じた式とのずれ（δ）は丸めずにそのまま出す（`run -v`）。

## 🟦 断面図（溜まる・流れる・遅れる）

```bash
uv run python fluidfullspeedsimulation.py draw flume-hopper --frames 4 --rows 4   # じょうご → 管（遅れ 3 秒）→ 樋
uv run python fluidfullspeedsimulation.py draw flume --fluid namekian-water         # 🟩 Namekian Water（🏷️ 仮置き）
uv run python fluidfullspeedsimulation.py draw flume --fluid ethanol --fill 🟥      # 色は自由（--fill で何色でも）
uv run python fluidfullspeedsimulation.py draw moniac-rig --frames 4               # 🏛️ 試験台の槽が溜まっていく
uv run python fluidfullspeedsimulation.py draw somen-flume --cap-mode rim          # 🥢 竹の縁を縦いっぱいに
uv run python fluidfullspeedsimulation.py draw bamboo-wear                         # 🎋 節の手前にたまり、うどんが止まる
```

```
⏱️ 30.05
🍜⬜🟦🟦🟦🟦🟦🟦🟦🟦🟦⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜
🟦⬜🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦⬜⬜⬜⬜⬜⬜⬜⬜⬜⬜
🟦⬜🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦⬜⬜⬜⬜⬜⬜⬜⬜⬜
🟦⬜🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦🟦⬜⬜⬜⬜⬜⬜⬜⬜⬜
```

- 列: じょうご ｜ 管（🍜 が中にいる間 ⏳）｜ 樋の 20 区間。行: 高さ（上が深い）
- 📏 段数 = round(水深 / 縦いっぱい × rows)。1 段 = 縦いっぱい / rows なので、段数から水深を ±½ 段で読み戻せる（E09 で確かめる）
- 縦いっぱいは `--cap-mode run-max`（走らせた中の最大。樋は全区間で共通）か `rim`（樋の縁）。`rows` の既定は 8
- 🍜 は水面の段に置く。ブロックの色は自由（`--fill`）。既定は 🟦 水・🟨 エタノール・🟧 メタノール・🟩 Namekian Water・🟪 Poison Water

## 🟦 ドット絵の水

| 💧 dam-break | 💧 moniac-drain |
|---|---|
| ![dam-break](docs/images/dam-break_water_s0_t0800_r01.gif) | ![moniac-drain](docs/images/moniac-drain_water_s0_t0600_r01.gif) |

```bash
uv run python fluidfullspeedsimulation.py draw moniac-drain --fluid ethanol --scale 8 --out-dir build/pixel \
    --file-prefix demo- --file-suffix _try   # 🎞️ GIF + 🖼️ PNG + 🧾 受領記録 JSON
```

- 🧱 `pixelwater.py`: 1 行 = int 1 つのビット盤。落下（真下 → 斜め下）は毎ティック、横の広がりは流体ごとの整数スケジュール（base_rate × ν_water/ν、Bresenham 型の積算）。⚠️ 見た目のたとえで、粘性流れの物理モデルではない
- 📛 名前: `<prefix><場面>_<流体>_s<seed>_t<ティック>_r<NN><suffix>`。`r<NN>` は同じ条件の何回目か（既存なら次の番号）。隣の `<name>.json` に plan（再現の引数）と bytes・sha256
- 🎨 色は Ironmate の色 ID（💧 C02・🍶 C04・🧴 C01、壁 C07）と派生規則。`to_ironmate_artifact` で `ironmate-pixel-grid/1` に出せる
- 🎞️ GIF は標準ライブラリだけ（LZW・前フレームとの差分矩形）。`scale` は 1 マスのピクセル数（既定 8）

## 🧮 ビット計算

- 🔢 Q 形式（`FRAC_BITS = 48`）。掛け算は `(a*b) >> 48`（床）、床の偏りは n ステップで高々 n ulp
- ✂️ `take`: 取り分を切り出して端数は残りへ → 分けた合計が元と完全一致
- 🔀 1 次遅れ（水槽 1 段）は `k = 2^−s` なら `(x − y) >> s` のシフトだけ
- √ Torricelli・堰は `math.isqrt` の整数平方根、樋は「セル体積 → 流量」の表を二分探索で引く
- 🎲 試す値は A×10^−N の A と N を乱数で選ぶ（`fixedq.random_scales`、Aoi の `draw_coefficients` と同じ作り方）。10^−3 刻みの決め打ちにしない

## 🧪 流体（20 °C）

| | ρ [kg/m³] | μ [mPa·s] | ν [mm²/s] | 出典と幅 |
|---|---|---|---|---|
| 💧 water | 998.21 | 1.0016 | 1.003 | List of viscosities / chem-casts（1.0016–1.002） |
| 🍶 ethanol | 789.45 | 1.20 | 1.52 | Lange（密度）、infobox 1.2 / CRC 1.144（幅 1.144–1.20） |
| 🧴 methanol | 791.1 | 0.587 | 0.742 | Flottweg の表（ρ = μ/ν） |

| 🟢 namekian-water | 998.21 | 1.0016 | 1.003 | 🏷️ 仮置き（名前と 🟩 だけ。物性は水と同じ） |
| ☠️ poison-water | 998.21 | 1.0016 | 1.003 | 🏷️ 仮置き（名前と 🟪 だけ。物性は水と同じ） |

`fluids.FLUIDS` に 1 行足せば流体を増やせる（出典と幅も一緒に書く）。🏷️ 仮置きの流体は描画では選べるが、実験（E01〜E11）の流体一覧には入れない。

🧊 水温: `get_fluid(key, temp_c)`・`Network(..., temp_c=10)`・spec の `"temp_c"` で μ・ρ が温度表（`fluids.TEMP_TABLES`）から動く。温度での変わり方を表から取り、20 °C では上の表の値にそろえる（比を掛ける）ので、温度を指定しない計算は 1 ビットも変わらない。

| | 範囲 | 出典 |
|---|---|---|
| 💧 water（🏷️ 仮置きも同じ） | 0–40 °C | [Anton Paar Wiki](https://wiki.anton-paar.com/en/water/)（IAPWS 2008、2–40 °C）、0 °C は chem-casts |
| 🍶 ethanol | 0–40 °C | μ: chemeurope Ethanol data page（CHERIC）、ρ: Lange の直線。Kaye & Laby（0 °C 1.873・30 °C 0.983 mPa·s）とは 3〜9 % ずれる |
| 🧴 methanol | 0–25 °C | [Anton Paar Wiki](https://wiki.anton-paar.com/en/methanol/)（Kaye & Laby ほか） |

## ⚠️ 限界

- 入口・出口損失、温度変化、表面張力は入れていない。Cd（0.61・0.62・0.58）と Manning の n（0.010）は代表値
- 陽解法（同時更新）。樋は Courant ≤ 1 を確かめる。定常水位には流入をステップ頭に入れる分 δ ≈ dt/2τ が残る
- 管の遅れは初期水位の流速で決めた栓流。トレーサーは先入れ先出し（水槽では平均滞留時間に当たる）
- 層流の樋は v = gSR²/(βν)（β: 長方形 3・半円 2）の近似
- 🥢 継ぎ目のたまりは手前の 1 セルにまとめる（上流へのせき上げは解かない）。曲がりの持ち上がりは Δh = v²B/(gR) の近似
- 🧊 越流（継ぎ目・節）の式は粘度によらないので、冷やしても引っかかりは変わらない。実際は越流水深が mm 程度だと粘性・表面張力で Cd が下がり、冷水の方が引っかかりやすい可能性があるが入れていない
- 🎋 割れ目は Cd 0.61（薄い縁の穴）。竹の厚み（数 mm）を通る短い管とみれば 0.82 程度まで上がる（E16 の summary に両方）。mm の水深では表面張力で漏れが止まることもあるが入れていない
- 🍜 麺の太さは目安（機械製麺の区分: そうめん < 1.3 mm ≤ ひやむぎ < 1.7 mm）。引っかかりは「越流水深 < 太さ」だけで決める

## 🗺️ この先

- 🎨 図（文字を入れるならメイリオで統一）
- 🍜 流しそうめん: 麺の長さ・からまり、竹の向きを変える分岐、箸ですくう取り出し
- 🟩 Node 側の鏡合わせ実装（BigInt で同じ整数列になるか）

## 📚 出典

- [Phillips Machine — Wikipedia](https://en.wikipedia.org/wiki/Phillips_Machine)（1949、Lancaster 爆撃機の余剰部品、約 2% の精度、Financephalograph）
- [Ethanol (data page) — Wikipedia](https://en.wikipedia.org/wiki/Ethanol_(data_page))
- [List of viscosities — Wikipedia](https://en.wikipedia.org/wiki/List_of_viscosities)
- [Flottweg: Dynamic viscosity](https://flottweg.com/wiki/separation-technology/dynamic-viscosity)
- [chem-casts: Water viscosity table](https://chem-casts.com/data-tables/water-viscosity-data-table)
- [ふるなび: そうめんとひやむぎの違い（食品表示基準の太さ区分）](https://furunavi.jp/discovery/knowledge_food/202504-somen_hiyamugi/)
