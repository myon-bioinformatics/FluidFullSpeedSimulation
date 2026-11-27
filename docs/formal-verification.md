# 保存則の形式検証

FFSS の実行時依存は Python 標準ライブラリだけです。Lean は任意の検証用ツールであり、
通常の実行や Python テストからは呼び出しません。共通の実行・結果判定は親リポジトリ
[`myon-bioinformatics/lean_verify.py`](https://github.com/myon-bioinformatics/myon-bioinformatics/blob/6c2558c4251698c50bc07bad86fb49f902b2738f/lean_verify.py)
で管理し、FFSS は対象となる命題と Python 側の接続テストを管理します。

## 証明する範囲

`formal/Conservation.lean` は量を数学的な整数 `Int` として扱います。
整数の 1 は Q48 の 1 刻み (`2^-48`) に対応します。
固定小数点の移動量が決まった後の、在庫の増減に関する命題です。

| 命題 | 内容 | Python 側の接続先 |
| --- | --- | --- |
| `FFSS.transfer_conserves` | 移動元から引いた整数を移動先に足すと合計が変わらない | `Network.step` の遅延なしの移送 |
| `FFSS.transfer_nonnegative` | 移動量が非負で在庫以下、移動先も非負なら、移送後の両方が非負 | 在庫不足時の上限処理と複数出口の配分 |
| `FFSS.queue_push_conserves` | キューの新在庫と出力量の合計は、旧在庫と入力量の合計に等しい | `DelayLine.push` と `in_transit` |
| `FFSS.delayed_transfer_conserves` | 水槽・キュー・到着待ち量を含めると、遅延付きの移送でも合計が変わらない | `Network.step` の `line`・`pending` と `total_q` |
| `FFSS.delayed_transfer_nonnegative` | 移動元・キューの上限と到着量の非負条件を満たせば、更新後の各在庫が非負 | 遅延付きの移送に対する境界条件テスト |

遅延付きの命題では、`queued` はキュー内の量、`pending` は今のステップに到着する量、
`released` はキューから出て次の到着待ちになる量です。次の等式を証明します。

```text
(source - sent) + (destination + pending)
  + (queued + sent - released) + released
  = source + destination + queued + pending
```

`pending` も保存量に含める必要があります。遅延が 1 ステップならキューの長さは 0 でも、
到着待ち量は存在します。

これは **Python 実装全体の形式的な正しさを証明するものではありません**。
流量の物理モデル、浮動小数点からの変換、配分アルゴリズムの全入力での正しさ、
複雑な流入補正、成分モードなどは、このファイルの証明範囲外です。
非負性の定理には仮定があり、実装がすべての状況でその仮定を満たすことも別の課題です。

## 三つの検証を組み合わせる

1. **数値実験・回帰テスト**: 既存の `test_network.py`、`test_fixedq.py`、
   `test_experiments.py` が、実際の回路や実験結果を確認します。
2. **境界条件・反例探索**: `tests/test_formal_bridge.py` は小さい整数の状態を列挙し、
   空の水槽、在庫を超える要求、整数の端数、Q48 や浮動小数点の精度を超える大きな整数、
   キューの長さ 0、複数出口、循環と遅延を含む固定シードの回路を検査します。
   有限個の入力によるテストなので、全入力を証明したことにはなりません。
3. **形式証明**: Lean が、明示した仮定の下で上記の整数の命題を検証します。
   `#print axioms` で各命題の依存する公理も出力し、専用 CI が判定します。

接続テストは `Fraction` から厳密な Q48 値を作り、実際の `Network.step` と
`DelayLine.push` の更新結果を確認します。模擬エンジンの結果を検査するものではありません。
xprobe などで探索範囲を広げる場合も、これらの接続点と仮定を維持します。

## 実行方法

通常のテストには Lean は不要です。

```sh
python -S -m unittest discover -s tests -p test_formal_bridge.py -v
python -S -m unittest discover -s tests -v
```

Lean は `formal/lean-toolchain` にバージョンを固定しています。elan を利用する場合、
同ファイルが適用される `formal` ディレクトリで実行します。Mathlib は不要です。

```sh
cd formal
lean -DwarningAsError=true Conservation.lean
```

結果を共通形式で保存するには、FFSS のルートから親リポジトリのラッパーを実行します。
以下は親リポジトリを隣の `myon-bioinformatics` ディレクトリに配置した場合です。

```sh
python -S ../myon-bioinformatics/lean_verify.py \
  --cwd formal \
  --required \
  --json-output build/lean-verification.json \
  --theorem FFSS.transfer_conserves \
  --theorem FFSS.transfer_nonnegative \
  --theorem FFSS.queue_push_conserves \
  --theorem FFSS.delayed_transfer_conserves \
  --theorem FFSS.delayed_transfer_nonnegative \
  Conservation.lean
```

`--cwd` は証明ファイルの相対パスと Lean ツールチェーンの基準です。
結果ファイルの相対パスは、ラッパーを呼び出したディレクトリを基準にします。
ローカルで Lean が未導入でも結果だけを記録したい場合は `--required` を外します。
未導入時は未検証として記録され、証明に成功した扱いにはなりません。

専用の Lean CI では共通ワークフローをコミット SHA で固定して呼び出し、
`--required` 相当の判定を使います。証明エラー、未導入、タイムアウト、公理の監査失敗は
CI の失敗です。通常の Python CI は引き続き独立して動きます。

命題と実装の変更は一緒にレビューし、Lean の成功結果と Python のテスト結果を別々に
確認してください。保存則の証明だけで、数値精度や物理モデルの妥当性まで保証はできません。
