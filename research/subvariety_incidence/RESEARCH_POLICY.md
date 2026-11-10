# 完全等方Log点：非被約接錐と周囲障害

2026-11-10 / v16。再開：ルートRESEARCH_POLICY.md → 本方針 → NOTES.md → papers/isotropic_log_tangent_cone_ja.tex。

## 一般原則・保存許可

研究価値、正確性、理解の深化、新規性、実現可能性を総合して方針を変更する。証明付き・今回の検算・未証明・既知の入力・新規性未確認を区別する。日本語LuaLaTeXの現行原稿と必要な入力を保持し、ノートに証明全文・旧版・実況を重複させない。ユーザーの保存許可の現行運用はGitHub mainを唯一の正本とし、非force更新後に先端を確認する。PDF非追跡、Actions不使用、Libraryは実際のpush失敗時のみ緊急退避。既存Libraryファイルは無断削除しない。生成・検算・組版・保存を分けて報告する。

## 現在の中心定理

q=u0²+2u1w0+2zw1、i_(P0)Ω=dq∧i_ξ(du0∧du1∧dz)、ξ=(1,1,-2)。30次元接空間の基底は新原稿第2節とproducerに明記。

s_i=t_(20+i), a=t27, b=t28, c=t29として、**全次数二初期イデアル**は
I2=(s0*c-s2*a, s1*c, c*(s2+c), s3*c, s4*c+s6*a, s5*a, s5*c)。
この7次元性は全150一次方向・350係数・465組の二次項により厳密に確認した。h=a*(s2+c)は次数二関係に入らず、h²はI2に入る。そのため**真の接錐でもh≠0、h²=0**であり、接錐は非被約。

二次障害スキームの被約台は次元28/26/24の三成分。最大次元成分だけは実接錐でもV(a,c)で、一般長さ1。対数型源の独立な28次minorと合わせ局所次元28、Hilbert–Samuel重複度1。これらは147成分完全性を仮定しない。

固定Q,σ延長のac関係は全周囲では明示R_acにより二次まで消去可能。P0+aK1+cK3+acR_acはC[a,c]/((a,c)^3+(c²))上でPoisson。ただし純c方向のc²障害は消えず、全次数への持上げも未証明。

既存の25次元長鎖rank4成分へのPGL5同定、28次元Log1112成分との交叉、非等次元・非CM・depth≤25は同伴原稿v11に保持。これらの成分性入力と新しい独立接錐計算を区別する。Q3一般三型、共鳴境界、cubic・del Pezzo剛性、P3対称性、147分類の必要数学資産は削除しない。

**新しい構造定理（v13）**：唯一の28次元rank2 Log1112成分はpで滑らか。e(R)=1の最大次元成分への結合公式でその成分自身の重複度が1となり、excellentな局所整域にNagata判定を適用。T_p Zlog=V(a,c)、形式局所環は28変数冪級数環。20次元のPGL5軌道はrank2とrank4の両成分に含まれ、交叉局所次元は23以上24以下。コード check_isotropic_orbit20.py と証明書 certificates/isotropic_orbit20.json は全25個のgl5微分を独立検算する。このsmoothnessは概念的帰結であり機械計算だけから推論していない。数値下界23は局所30変数の一般交叉次元不等式からも従うが、q_vとxi_tとb核方向を使った明示的な3パラメータPoisson共通族を構成した。rank2はLog1112のepsilon留数極限、rank4は非定数g長鎖源とアフィン線形座標変換で全六括弧が一致。gl5と三追加方向のrank24、非零小行列式10871635968、Eulerゲージと合わせたrank48が源像の23次元を保証する。原稿thm:log-long-intersection23、独立producer check_isotropic_intersection_dim23.py。交叉の正確な次元23/24、scheme構造は未解決。


**新しい普遍延長（v15）**：A=u1∂u0−u0∂w0の下で滑らかなq=u0²+2u1w0+2zw1+αu1²+βu1z+γu1w1+δz²+ηw1²を動かしても、P(q,ξ_t)+λK_A(q)は全パラメータでstrict Poisson、Pfaffianは(3λ/2)u1²z dq。λ=0のLog1112対数留数表示とλ≠0の長鎖非定数gのアフィン線形正規化も全係数で確認。さらにw0'=w0+(αu1+βz+γw1)/2、w1'=s(w1+κz)、ηκ²−2κ+δ=0、s=1−ηκ≠0によりqは一変数q_η'（η'=η/(1−δη)）に帰着する。したがって固定Aの五係数を探索し直す価値は低く、24次元交叉方向があれば固定Aを破る必要がある。証明papers/isotropic_nilpotent_pencil.tex、独立検算check_universal_nilpotent_pencil.pyとcertificates/universal_nilpotent_pencil.json。これは分類済みの長鎖成分性を入力にし、元の局所環の被約性は解決しない。

**固定冪零場の二軌道と収縮（v16）**：q_η=u0²+2u1w0+2zw1+ηw1²、P_(η,t)=P(q_η,ξ_t)に対し、g_ρ=diag(1,1,ρ,1,ρ^-1)はAと可換、det=1、射影等式g_ρ·[P_(η,t)+λK_η]=[P_(ρ²η,t)+ρλK_(ρ²η)]を全10係数で証明。ρ→0で完全等方rank2点[P_(0,t)]に収縮。基底三形式の共通核Wへの制限はrank(q|W)=0（η=0）または1（η≠0）で、Aと基底三形式の直線を保つqは二軌道に分かれる。λ≠0のη/λ²はこのC*群の不変量だが全PGL完全不変量とは未確認。証明thm:nilpotent-torus-attraction、独立producer check_nilpotent_torus_attraction.py、certificates/nilpotent_torus_attraction.json。

## 次の核心・着手方法

最優先は固定Aと三形式の直線を保つ二軌道の外の変形による交叉次元23/24の判定と、全初期イデアルの三次以上。線形化の95 normal方向を二次まで解き、30 tangent方向上の三次残差とI2のsyzygyを計算する。二次模型の26・24次元成分は真の成分と未同定であり、新しいPoisson成分として数えない。
次にhの実接錐上の台・annihilator、既知25次元成分の接錐とscheme交叉を調べる。元の局所環Rの被約性と全成分の分類は未解決。接錐非被約からR非被約を推論しない。局所分類の価値に応じてより一般的な障害機構やdel Pezzo例外へ方針を変更する。文献上の優先権は未確認。

## 再現

新原稿 papers/isotropic_log_tangent_cone_ja.tex は独立6頁。コード code/subvariety_incidence/check_isotropic_quadratic_cone.py と check_isotropic_intersection_dim23.py、証明書 certificates/isotropic_quadratic_cone.json と certificates/isotropic_intersection_dim23.json。新しい手動入口 bash scripts/verify_github_source.sh tangent-cone。SymPy 1.14.0の全有理計算と6頁LuaLaTeX二回組版を実行済み。旧11証明書・147分類の全面再実行は今回行っていない。
