# nonsquare lazy quartic oracle の独立 source・DAG・係数検証

2026-11-03。元の134式と G6 を直接読み、βを代入して t で二次へ斉次化し、
q_fixed=1 として全135個の J、B を新たに作る。著者の source 抽出関数や
producer を実行しない。保存された全 sparse J/B と一致し、JT=0、選択
34行の JU=I、全3484個の core atom の分母が既知 aFBC のみと検査する。

全 source J/B の局所重み、T/U の重み、nuisance 二座標の恒等零も検証する。
全1373 node の順序、重み、ordinary degree を調べ、元の J/B と正確な core
算術から w2、w3、w4、各 residual の DAG を再構築して同じ node に到達する。
core の全パラメータ四次行列を展開して検査したのではない。

選択式の重みは1〜4、normal の重みは全て正である。従って5次以上の
ordinary monomial は選択式の重みに入らず、w2+w3+w4 はその全34式の厳密な
多項式 graph である。残る高重みの式を零とする主張ではない。選択外の
元の weight4 行33本とその root の番号を全て source に照合する。

元の行77、131の H² だけを別に検査する。全 light を零、normal v5=H とし、
元 source から w2_H² を直接解いて `B(H,H)+J*U*w2_H²` を作る。一方で DAG
を独立の H 一変数係数配列で再帰評価する。この二つと、旧凍結 source の
二つの有理係数公式が QQ(a,e,d) 上で厳密に一致する。ordinary3/4 成分には
H² は入らないことも確認する。最後に公開 coefficient API の二つの query
が同じ値になることを照合する。新しい core 点や有限体は使用しない。

元のシリアライズ済み DAG は1373 node、二 query 後のメモリ内は1407 node。
旧 API receipt の original count が共有 list の変更後を指すという凍結済み
clarification は正しい。ファイル自体のハッシュは query 前後で不変である。

以上の範囲で oracle と二係数 API を受理する。全122×64行列の rank64、
quartic 全消去、ambient 除去は依然この検証の結論ではない。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/s61_nonsquare_lazy_quartic_oracle_review_20261103_v63.py --output work/nonsquare-oracle-independent-review
```

60秒・2.2GB上限。元の producer は再実行しない。
