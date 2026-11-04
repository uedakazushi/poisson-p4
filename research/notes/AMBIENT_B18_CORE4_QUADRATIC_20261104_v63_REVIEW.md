# CORE4二次initialとsource三次補正消去の独立監査

2026-11-04。対象は `CORE4_QUADRATIC_PREFLIGHT` と `CORE4_QUADRATIC_INITIAL` の凍結二packetである。actual N-pair hypersurface、CORE4 reduction、共通physical R0二次模型の既受理独立gateを使用する。source solver、source高次jet、著者worker、GB、追加点は使わない。

## 元sourceと四変数二次近似

typed-core4のscalar DAGが元actual-sourceの全prefixをそのまま保持することを確認する。二族48式ずつと二つのmonic式がactual98系であること、未知数が元96係数とx27,x61であることをliteralに結ぶ。独立自由多項式演算で全十二output差が同じphysical aA+bA tauA−aB−bB tauBであることも確認する。二族のtauは同一視しない。

保存された両57 source二次cacheを、各族F9のcommon R0逆写像、二つのmonic Schur二次写像、四変数(x112,x117,tauA,tauB)へのtangent制限の順に、別のQQ dictionary多項式合成で戻す。全107 source/parameterの15係数jetを照合する。新しいsource解を求めたのではない。

七つの修正F'は物理rows67,72,63,64,68,73,74で、row72だけrow60を引く。この七式が**full九変数Schur環で**二次以下零であることを保存E4二次模型から確認する。四変数に制限した後だけの零性を全体へ外挿しない。

## lambdaによる未知三次係数の消去

必要な元DAG祖先17085個について、原点値と98未知数への疎な一次微分を独立に計算する。全98×98 Jacobian J、七行B=dF'/du、保存lambdaをQQで照合し、全686項目のlambda*J=Bを確認する。unit determinantは−12524124635136/78125である。

この例では保存lambdaの七行もBの七行も実際に全て零であることを別に確認する。従ってsource三次補正が寄与しないことは直接にも分かる。一般の式として、u=u_[2]+u3+O4、E(u,c)=0ならJ*u3=−E3(u_[2])であり、真のF'3はF'3(u_[2])−lambda*E3(u_[2])となる。未知u3を真のsource解で零と仮定したものではない。

## 全245係数とrank2

四変数次数≤3の35単項式をkeyとする独立QQ dictionary interpreterで、17085祖先を再評価する。全98式の二次以下零と保存三次残差、七raw/corrected三次jetの全係数を確認する。unit逆元は有限幾何級数を用い、著者の固定畳込み表は使用しない。

L1=5*x112−20*x117+24*tauA+6*tauBと、保存された七二次式の積を作り、七×35=245係数すべてで補正済み三次式へ戻ることを確認する。二次係数行列の独立QQ rankは2で、最初の二行のみ非零、残る五行は零である。選択2×2 minorは703/82944で非零となる。

## 二次initialをactual coreへ移す論理

既受理の全次数F'=L*G'と、full Schur環でのF'二次以下零から、G'の一次以下は零である。従ってF'3=L1*G'2となり、L2/L3は不要である。五つの選択Gをunit IFTで消すと、消える五physical座標は四core変数に対してorder≥2である。G'に一次項がないため、その二次以上の加速度はG'2を変えない。row72−row60という修正も、選択G60=0上では元のcore関数と一致する。

よって確認した七二次式は真の七core関数の二次initialである。七関数がorder≥2なので、そのidealの二次部分全体はこのrank2の定数係数線形包となる。高次関数や相殺からさらに高次initialが生じる可能性は残る。

これは全core ideal、socle、nilradical、E4またはD0のdepthの計算ではない。零だった五二次式が全次数で零だとも言わない。source三次求解、N-pairの新しい高次係数、source四次への自動引上げは行っていない。全input/evidence hashと実測値は同名REVIEW/report.jsonに記録する。
