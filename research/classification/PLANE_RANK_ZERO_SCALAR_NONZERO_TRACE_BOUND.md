# scalar first symbol の全 T非零源：projective 像は二十三次元以下

2026-11-07。主義務 G、producer `plane_collision_review`。
**命題。** 実既約被約 quintic D を持つ rank 四 Poisson tensor が plane P 上零となり、
traceless co-Higgs first symbol が恒等零、scalar trace field T が非零なら、
その全 projective Poisson 像は次元≤23。
位数二・ordinary double・semisimple性・generic eigenvalueを仮定しない。
T0K0及び traceless非零でtraceとpointwise dependentなglobal-V源は本命題の外に保持する。

全 source と actual inverse normalization、正確な source 次元は
`PLANE_RANK_ZERO_SCALAR_FIRST_SYMBOL.md` で証明した。
本稿は三つの有限 exact operator 証明書をその全十五modelへ戻す。
独立監査前の producer 証明であり、新成分又は E帰属は主張しない。

## 1. 全 first-symbol 完全性と dimension の戻し方

全重複度 conormal 可換性から full35/30/9 source を使用する。
scalar 条件の全 weight零源は
\(\Pi_0=T\wedge E_n+q_S\partial_z\wedge\partial_w\)、
\((T,S)\in\mathfrak{sl}_3\times\mathfrak{sl}_2\)、affine11である。
全 symbol の Jacobi/Pfaffian は零だが、それだけから higher-source の空性は推論しない。
normal conjugation、base conjugation、tensor scalarと全逆によって
10固定・4一parameter・1二parameterが全 T非零を覆う。
全30/9 unknown は各 normalization で保持される。

各固定 first symbol に対し全350係数から D1,D2を作り、
\(d=\dim\ker D_1+\dim\ker D_2\) を使う。
任意 P1 上で P2 は同一 D2 の inhomogeneous fiberなので、全 higher Poisson fiber は d以下。
affine symbol source dimension \(d_T+d_S\) を戻すと全像は

\[
 d_T+d_S+6-1+d.
\]

finite parameter wallではその条件のcodimensionを sourceに戻す。
fixed normal-form kernel をそのまま ambient source dimensionと混同していない。

## 2. 固定十型：九つの rank bound と一つの全源因子

契約 `jobs/plane_rank_zero_scalar_fixed_weights.json`、不変 report
`reports/CLASSIFICATION_PLANE_RANK_ZERO_SCALAR_FIXED_WEIGHTS_20261107.json`。
1.211947秒、全算術 PASS。
全150 P0、全350 Jacobi/Pf0、全350×30・350×9、完全kernel、nonzero minorを保存した。

| fixed T/S | dimker D1 | dimker D2 | 全像linear upper |
|---|---:|---:|---:|
| N3 / 0 | 8 | 3 | 22 |
| N3 / E | 6 | 0 | 19 |
| N3 / H | 6 | 0 | 20 |
| N2 / 0 | 16 | 6 | 31（この bound 単独では未閉鎖） |
| N2 / E | 7 | 0 | 18 |
| N2 / H | 6 | 0 | 18 |
| D / 0 | 6 | 0 | 16 |
| D / E | 6 | 0 | 18 |
| D+N2 / 0 | 6 | 0 | 18 |
| D+N2 / E | 6 | 0 | 20 |

N2/0では T=b∂a、P0は∂aを共通因子に持つ。
全D1 kernelを読むと、∂aを含まない方向は∂b∧∂cだけで、normal-mixedの非∂a方向は零。
全P2もbase pairsだけなので、非線形higher equationを解く前から全tensorは

\[
 \Pi=\partial_a\wedge V+F_2\partial_b\wedge\partial_c
\]

となる。actual PfaffianはF2の倍数。
F2零ならrank≤2、F2非零なら二次×三次で、実既約被約quinticに矛盾する。
従ってこの一源は空であり、固定十源全てが≤22となる。
初回reportの `all_ten_fixed_G24_closed=false` は linear判定だけの結果として、そのまま保持する。

## 3. 全四一parameter源と全根

契約 `jobs/plane_rank_zero_scalar_one_parameter.json`、不変 report
`reports/CLASSIFICATION_PLANE_RANK_ZERO_SCALAR_ONE_PARAMETER_20261107.json`。
2.494044秒、全算術 PASS。
全 QQ(τ) kernel/minorと全 irreducible minor rootをexactに保存した。
今回の全factorはQQ上一次であり、数値近似を使わない。

| 全 source | domain | generic d / image | 全domain内wall最大image |
|---|---|---|---:|
| T=τD,S=H | τ≠0 | 6 / 19 | 19 |
| T=τD+N2,S=H | τ≠0 | 6 / 21 | 20 |
| T=diag(1,τ,−1−τ),S0 | τ≠1,−2,−1/2 | 6 / 19 | 23 |
| 同じT,S=E | 同じdomain | 6 / 21 | 20 |

前三角固有値の scale parameterは全 source 次元へ戻した。
前二源の全 selected-minor rootsは \(0,1/3,2/3,-2/3,-4/3\)。
τ0はそのJordan domain外であり、T0又はN2として別に保持した。
後二源の全 rootsは \(0,-1,1,-2,-1/2\)。
最後の三 root はrepeated semisimple源へ復元される。
S0のτ0,−1ではd11、fixedparameter source dimension七を戻して像23。
他の全domain内rootは表の弱いboundで十分。

outside-domain arithmeticを捨てず reportに保持し、現在の regular-Jordan theoremへ誤って含めていない。
T0は本命題そのものの外なので、ここで別作者の未監査T0証明を入力にはしない。

## 4. 全二parameter regularSS/H：十二lineと全五十七root

契約 `jobs/plane_rank_zero_scalar_two_parameter.json`、不変 report
`reports/CLASSIFICATION_PLANE_RANK_ZERO_SCALAR_TWO_PARAMETER_20261107.json`。
5.882620秒、RSS59,924KiB、全算術 PASS。
\(T=\operatorname{diag}(\lambda,\mu,-\lambda-\mu)\)、S=Hの全源を扱う。
regular域は \(\delta=(\lambda-\mu)(2\lambda+\mu)(\lambda+2\mu)\ne0\)。

global QQ(λ,μ) certificate はdimker D1=6、D2 injective、affine symbol11から全像22。
selected-maximal-minor productとδの全factorは次の十二本である。

\[
\begin{gathered}
 \lambda,\ \mu,\ \lambda+\mu,\quad
 3\lambda+2,\ 3\lambda+4,\ 3\mu+2,\ 3\mu+4,\\
 3\lambda+3\mu-2,\ 3\lambda+3\mu-4,\quad
 \lambda-\mu,\ \lambda+2\mu,\ 2\lambda+\mu.
\end{gathered}
\]

最後の三本はentirely δ0なので、repeated-Jordanの全一parameter源へ送る。
originのT0は本命題の外として保持する。
他九本の**全line**をexact affine parameter τで記述し、
full350行をQQ(τ)上で再計算して各lineのcomplete kernel/minorを得た。
全line genericでd6、symbolbase10を戻して像21。
各line自身のminorの全rootをさらにexactに検査し、合計57点証明書を保存した。
δ非零の全pointではd6、symbolbase9を戻して像20以下。
δ零のpointはoutside-domainとして正確に判定し、対応する別Jordan源へ送る。

global selectedminor非零のpatch、その全line、各line minor非零のpatchと全rootsでregular C2全体を覆う。
line同士のintersectionも、各line上のgeneric certificate又はそのexact rootに必ず入り、未検査の点を残さない。
finite pointを無根拠にgenericと呼ばず、nonlinear factorが現れれば停止する契約であり、今回は全factorが一次だった。
従ってこの全二parameter源は≤22。

## 5. immutable execution と最終結論

実行器は一つの reusable `code/classification_plane_rank_zero_scalar_weights.py`。
fixed十型→一parameter→二parameterはisolated modeとして加え、
各旧reportに埋め込まれたproducer sourceの全FunctionDef ASTが不変とruntimeで確認された。
旧report bytesは変更せず、旧jobのrunner bindingだけexplicitなcompatible-extension理由とhistoryを付けてrefreshした。
旧high-m/κpacket、shared fullrank runner、global/main/Gitは変更していない。

固定十源≤22、全四一parameter源≤23、全二parameter源≤22を合わせ、
**traceless co-Higgs 零・T非零の全scalar first-symbol源は≤23**。
高weight非線形Jacobiはこの必要sourceを増やさない。
T0K非零の別producer全source手証明、T0K0、traceless非零のdependent global-V源を混同せず保持する。
全G/Eの完了や個別成分帰属は本稿から主張しない。
