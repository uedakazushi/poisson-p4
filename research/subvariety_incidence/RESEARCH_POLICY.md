# Poisson部分多様体incidence：現行研究方針

2026-11-11 / v17。再開順：ルートRESEARCH_POLICY.md → 本方針 → NOTES.md → papers/poisson_subvariety_incidence_ja.tex → 関連共有TeX。

## 一般原則と保存許可

研究価値・正確性・理解の深化・新規性・実現可能性を総合し、現行の証明と課題を優先する。証明済み／計算確認／条件付き／予想／新規性未確認を区別する。詳細な証明は現行日本語LuaLaTeX原稿に保持し、研究ノートに全文を重複しない。GitHub mainが唯一の正本、非force pushと先端読戻し、GitHub Actions不使用、生成PDF非追跡、Libraryは実push失敗時の緊急退避のみ。旧版・履歴・中間生成物は保存しないが数学的内容や検証可能性は削らない。全面先行研究調査と外部査読は未完了。

## 上位目標と第一の課題

**上位目標**：Poisson不変部分多様体の幾何学的incidenceからP4の147被約Poisson成分の源・退化・交叉を構造的に説明する。

**第一課題**：滑らかな四次del Pezzo S=V(q0,q1)上の非零Poisson切断σ_hの全ambient延長を特殊な組まで調べる。一般の(S,σ_h)の全延長がrank2の被約一点であることは証明済み。全候補は P=J(q0,q1,h)/4+q0 A+q1 B、A,Bが定数反対称行列の20変数二次Jacobiスキーム。特殊な組でrank4が存在するか、普遍剛性か例外分類かは未解決。

## 証明された幾何的な検査基準

- **v17の必要条件**：C=S∩V(h)は(2,2,1)の反標準曲線、次数4、算術種数1、一般には滑らかな楕円曲線。Sを保存するambient Poisson πはC上で零となり、Pfaffian五次式 Dπ∈I_S I_C⊂I_C²。rank4ならその退化五次超曲面はCに沿って特異。H0(I_S(5))は65次元、H0(I_S I_C(5))は41次元で、C上の法線一次ジェットに24本の独立条件。楕円四次Cの族は20次元、固定Cを超平面切断に持つ滑らかなSのファイバーはA10の非空開。普遍的な証明は papers/del_pezzo_pfaffian_incidence.tex、独立数値検算は code/subvariety_incidence/check_dp_pfaffian_conormal.py および certificates/dp_pfaffian_conormal.json。**これはrank4延長の十分条件ではなく、双ベクトル全20変数のJacobi方程式に追加の24本の独立な式を課すものでもない。**
- **従来の結果**：smooth cubicは既知intrinsic分類から全点rank2一意延長（36次元Rat113 incidence）、smooth degree≥4 hypersurfaceの非零ambient不変性なし、一般del Pezzo被約一点と28次元Rat122、smooth quadric一般Rat/Log/Affの一点/A1/二点と特殊Log境界、単純hyperplaneの横断対称性。
- **分類との関係**：141 rank4・6 rank2という既存147成分リストを照合対象とするが、構造定理の証明で完全分類を循環的に仮定しない。文献上の完全な優先権と外部査読は未確認。

## 次に行う具体的研究

1. 反標準楕円四次零曲線Cを先に固定し、長さ10の滑らかなdel Pezzo liftとrank4 Poisson源の不変性を比較。147源の実PfaffianだけでなくPoisson**零スキーム**を照合する。一般元での不成立は特殊元の排除にならない。
2. 対角形 q0=Σx_i²、q1=Σλ_i x_i²（λ_i相異なる）で特殊hの全20変数Jacobiイデアルを検査。特別な線形化minorと実rank4解を区別する。観察にとどまる標本は証明済みと称しない。
3. 次にspecial quadric Q3の制限形式・直交性・strict stabilizer・modular場・全Jacobiファイバーを層別し、147源の説明へつなぐ。
4. Fano hypersurface、hyperplaneとdel Pezzoのincidenceを統合する幾何学的構造定理を目指す。一般次元FanoとP5の(2,2)Fano threefoldは未解決の長期目標。

## 保持する副成果と再現

完全等方Log点の非被約接錐・局所次元28・23/24交叉等の有効な成果は papers/isotropic_log_tangent_cone_ja.tex に保持し、主課題の実行を妨げる単なる三次計算の継続は優先しない。数値証明書は厳密な有限計算、普遍論証は本文。手動再現は bash scripts/verify_github_source.sh incidence。組版、計算、GitHubへの保存は別状態として報告する。一次文献ではLoray–Pereira–Touzet、Mayanskiy、Pym surveyおよびNasu 2026を適用範囲付きで比較する。
