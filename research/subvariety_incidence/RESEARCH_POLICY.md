# 完全等方Log点：非被約接錐と周囲障害

2026-11-10 / v13。再開：ルートRESEARCH_POLICY.md → 本方針 → NOTES.md → papers/isotropic_log_tangent_cone_ja.tex。

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

**新しい構造定理（v13）**：唯一の28次元rank2 Log1112成分はpで滑らか。e(R)=1の最大次元成分への結合公式でその成分自身の重複度が1となり、excellentな局所整域にNagata判定を適用。T_p Zlog=V(a,c)、形式局所環は28変数冪級数環。20次元のPGL5軌道はrank2とrank4の両成分に含まれ、交叉局所次元は20以上24以下。コード check_isotropic_orbit20.py と証明書 certificates/isotropic_orbit20.json は全25個のgl5微分を独立検算する。このsmoothnessは概念的帰結であり機械計算だけから推論していない。

## 次の核心・着手方法

最優先は交叉次元の厳密値と全初期イデアルの三次以上。線形化の95 normal方向を二次まで解き、30 tangent方向上の三次残差とI2のsyzygyを計算する。二次模型の26・24次元成分は真の成分と未同定であり、新しいPoisson成分として数えない。
次にhの実接錐上の台・annihilator、既知25次元成分の接錐とscheme交叉を調べる。元の局所環Rの被約性と全成分の分類は未解決。接錐非被約からR非被約を推論しない。局所分類の価値に応じてより一般的な障害機構やdel Pezzo例外へ方針を変更する。文献上の優先権は未確認。

## 再現

新原稿 papers/isotropic_log_tangent_cone_ja.tex は独立6頁。コード code/subvariety_incidence/check_isotropic_quadratic_cone.py、証明書 certificates/isotropic_quadratic_cone.json。新しい手動入口 bash scripts/verify_github_source.sh tangent-cone。SymPy 1.14.0の全有理計算と6頁LuaLaTeX二回組版を実行済み。旧11証明書・147分類の全面再実行は今回行っていない。
