# 型7七全源の帰属：三MPS源、既存quintic源、新しい23軌道成分候補

2026-11-07、PRIMARY E。全parameterの被覆は七全源定理、各源の全PGL飽和閉包が
指定anchorの22軌道閉包に等しいことは `HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_ORBIT_MODELS.md`
による。本稿は残る第一・第二・第五・第六・第七源をliteralな実Poisson曲線で認識する。
第三・第四源の22次元成分性は同orbit-model稿で扱った。
本producerの計算は全PASS。独立再現と新第五成分の全prior99比較は確認中であり、countは変更しない。

## 実曲線から全源閉包への橋

各源のcanonical点P0から、係数が有理数の多項式二ベクトルP(t)を明示した。
証明書には全10括弧を保存し、全350 projective Jacobi係数がQQ[t]上で恒等零、
P(0)=P0であることを直接確認した。発見に使った二法線のIFTや切片GBはこの証明の入力でない。

anchorの座標C* weightをσ、曲線parameter weightをμとすると、全monomialで

    diag(λ^σ)_* P(t)=λ P(t/λ^μ)

が成り立つ。μは第一・第五・第六源で1、第二・第七源で2。
複素数体では全t≠0が同じ実PGL/scalar軌道に属する。
従って、その軌道を含むPGL安定閉集合CはP0も含み、
全源定理により当該の全parameter源のPGL飽和閉包を含む。
一点だけの帰属から任意parameterへ飛躍するのでなく、既済の全源22軌道閉包定理をここで用いる。

## 第一・第二・第七源は既知MPS支持に入る

三曲線ではh-Casimirを保持する。h-Casimir代表の斉次発散Zの行列Nを全係数で計算し、
特性多項式を順に

    第一源： u(t−u)(t+u)(2t−u)(2t+u),
    第二源： u(t−u²)(4t−u²),
    第七源： u(4t+u²)(16t+u²)

と得た。t≠0では全てregular semisimpleで、非零固有値比は±1,±2、
公開semisimple入口表の型12である。
実quartic Pfaffianのt=0は元の全源の既約quarticと正確に一致する。
可約quarticの集合は閉なので、一般tでも既約である。
従って一般曲線点は既済 `HYPERPLANE_QUARTIC_SEMISIMPLE_SUPPORT.md` の全型12支持定理で
有限個の既知MPS成分の和集合に入る。その和集合の閉性とPGL安定性により、
三つのanchor、さらに第一・第二・第七全源の全閉包が被覆される。
一意MPS成分や特定T番号までは主張する必要がない。

## 第五・第六曲線の非零点は滑らかな23軌道成分を与える

全非零tが同一軌道なので、構造的な代表t=1についてのみ計算した。
全150方向/全350行D、実GL5方向・Euler鉛直・scalarを合わせたMに対し、両点とも

    rank D=102,  rank M=48,  D M=0

となる。非零102/48小行列を保存した。
射影接次元と実軌道次元はいずれも23であるから、各点は滑らかで成分は唯一、
その23軌道閉包が全Poisson空間の既約成分である。
P0を経由する実曲線と先の全源22閉包定理により、第五・第六全源が対応する成分に含まれる。
任意標本点の探索は行っていない。

## 第六全源は既存99番目のquintic成分に入る

Euler代表を `P−(E∧div P)/5` にして発散を零にし、完全なprojective stabilizerを求めた。
その固有値は比例を除いて `(−2,−1,0,1,2)`、全て単純である。
証明書の有理可逆固有座標行列Bを使うと、primitiveな九括弧モデルが得られる。
さらに `D0=diag(1,1,−1,1,−1/48)` で押し出すと、
`CLASSIFICATION_GA9_DISTINCT_COMPONENT_AUDIT_20261107.json` の九括弧一般元の−1/12倍へ
全150係数で正確に一致する。元のprimitive点からの行列はD0 B^{-1}、scalarは−125/2である。
従って第六曲線の全非零点は既存
`R4_GA9_DISTINCT_IRREDUCIBLE_QUINTIC_23` の密軌道に属し、
第六全源はその閉包に含まれる。

## 第五の簡約一般元と複素絶対既約性

第五点にも同じ全周囲階数証明がある。発散零代表のprojective stabilizerは完全計算で一次元、
半単純で固有値の重複度は `(2,2,1)`、固有値は `(-512,-512,128,128,768)` と取れる。
有理可逆固有座標への変更後、primitiveな一般元は次の九括弧である。
座標は `(h,x,y,z,w)`、省略した{h,x}は零、反対称性を補う。

    {h,y}=−6h²,         {h,z}=3hx,
    {h,w}=−3hy+162hz−45xy,
    {x,y}=9hx,          {x,z}=−2x²,
    {x,w}=−5hz+2xy−243xz,
    {y,z}=−9hz+3xy,
    {y,w}=6hw−3y²+567yz,
    {z,w}=2xw+7yz−243z².

実五次因子は75Fで、

    F=h³z²+h²wx²+h²xyz+27hx²yz+3x³y².

Fはwの一次式で先頭係数h²x²を持つ。
残項Rの制限は `R|h=0=3x³y²`、`R|x=0=h³z²` なので、
複素多項式環C[h,x,y,z]でも係数h²x²とRのgcdは1である。
Gaussの補題からFはC上で既約。
QQ上の因数分解だけから絶対既約を推測していない。

`∂w F=h²x²` と全微分への直接代入から、被約特異支持は

    {h=x=0} ∪ {h=y=0} ∪ {x=z=0}

という三つの射影平面である。
実際、h=0,x≠0では∂y F=6x³yからy=0、x=0,h≠0では∂z F=2h³zからz=0となる。
逆包含は全五偏微分が各平面上で零であることから従う。
既約三次元超曲面に二次元特異支持があるのでFは非正規である。

この23成分は既存99番目とは密軌道のprojective stabilizer固有値重複度で区別できる。
既存99は全単純、第五は2,2,1であり、半単純Lie生成元の非零倍・共役でも変わらない。
別の区別として被約特異平面数も既存99の二面に対して三面である。
他の既存98項目との比較には、accepted Ga9 distinct auditの全同次元成分が閉可約五次locusに入る
という全一般源の証明を使える。異次元の既存全成分とは成分次元で区別する。
この全prior99比較と独立算術レビューの受理・台帳更新は統合担当が行う。

## 保存した証明書

`CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_CURVE_RECOGNITION_20261107.json`、
runnerの `--nil22-long-quartic-curve-recognition`、
`jobs/hyperplane_quartic_nil22_long_quartic_curve_recognition.json` を使用する。
literal曲線の係数はrunnerに固定し、公開orbit-model/七全源証明書と既存semisimple/99証拠だけを読む。
private IFT、jet、GB、workファイルは実行時入力でない。
旧37関数はAST不変であり、旧producer bytesもそのまま保存される。
