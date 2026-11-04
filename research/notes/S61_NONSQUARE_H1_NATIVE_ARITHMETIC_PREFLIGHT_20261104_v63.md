# H1四次係数の演算方式preflight（実時間profileではない）

113行cacheでの単元検査後、最初のquartic係数評価が遅かった原因候補を
既存実装のsourceから限定する。SymPy 1.14.0のQQはこのvenvではPythonMPQで、
python-flint moduleは無い。既存activate済みnative FLINT 3.1.0は利用できる。
環境の再セットアップやpackage導入は行わない。

`FracElement.__mul__`は分子分母を掛けて`new`し、`__add__`も共通分母が同じ
場合を含めて`new`する。`new`は毎回`PolyElement.cancel`を呼ぶ。QQ係数では
分子分母の有理数分母を払ってZZ上の多変数cofactor/GCDを計算し、QQへ戻す。
従って各DAG nodeで一般の多変数GCDが発生する経路がある。これはsourceの
callchainの観察であり、CPU時間の大半がGCDだったとの実測主張ではない。

入力2378 atomには228種類の分母があり、非monomialは52種類、最大174項・
全次数18・有理数係数1172bit。scalar nodeは174022個である。次の有限gateは
この228分母だけを既知a,d,F,B,C,Jへ正確に分解し、再構成を確認すること。
全て通れば、native polynomial numeratorと既知単元の指数vectorによる演算で
一般の多変数GCDを避けられる。加算では単元指数の最小値へ揃え、乗算では
指数を足し、必要な既知因子だけexact divisionで除く。最後のcross-product
恒等式は元scalar回路との一致を証明する。228分母分解や四次係数再評価は
このpreflightではまだ実行していない。
