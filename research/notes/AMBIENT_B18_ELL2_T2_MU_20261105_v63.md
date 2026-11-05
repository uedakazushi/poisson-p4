# honest T2 の実最小生成数

入力は別packet `ELL2_T2_MU_INPUT_GATE` のみで閉じる。元の完全132対1186表示、
27個のX0単位列、G3(0)=0、全102定数に単射なG2linearを用いる。
741個の定数zero-origin組合せごとに、実Xeff1の九方向から
G2linear Y0=Mnew0 Xeff1 の全234係数を確認する。未知M1はXeff0=0で消える。
旧八方向は独立検証済みの配列を再利用し、新ell方向だけ保存FULL9係数のds/5から
合成する。nativeはS0*cNative=0を毎回再確認してからL0*S1ellを使い、N1を作らない。

元のY列を一つずつ処理し、独立列、固定pivot、正規化basisの元列表示と全741列の
factorを保存する。全factorを元Y列に戻して照合する。固定元minorの行交換符号と
非零QQpivot積がrankの下側証明、全741factorが上側証明である。
最小生成数は正確に105-rank(Y0)。全長17..20やKのtorsion指数はこの計算だけでは
判定しない。新しいphysical/implicit/gauge/N1係数は計算しない。

上限は60秒、2GB、常駐展開QQ150000、384bit、RAW32MiB。各完了列をphaseへ記録する。
gateの最初のfree4読取り漏れは別directoryへ保存され、値workerの入力にしない。
