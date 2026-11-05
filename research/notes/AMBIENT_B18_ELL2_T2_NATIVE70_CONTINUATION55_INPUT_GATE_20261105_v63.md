# 未計算55例と全70 lower rank の固定 source gate

最初の15例は immutable な実値 prefix と独立 witness563 peer から読み、再計算しない。
今後の応答再生は case619 index564..618 の55例に限り、同じ C72・元 source・全 carry を使う。
entry の128 bit、temporary385 bit、numeric cap51263とRAW8305349は変更しない。
明示所有は元の49734に cached prefix二重510、rank scratch540、control100を加えた50884である。

rank は一般の有理 Gaussian 消去でなく、既知の非零 witness を必ず含む bordered minor を使う。
元の四行共通分母Dで整数化し、|v_i|≤U_i と既知 witness W_i を利用する。
サイズsの全 determinant の絶対値と途中和は (s−1)!Σ_i|W_i|Π_(j≠i)U_j で抑えられ、全15 row subset の最大は374 bitである。
これは385 bit temporary以内である。原有理係数と row/column index、非零 minor の literal recipe を保存する。
大きな評価済み minor や span係数を128 bitの committed データとして保存しない。
全70の固定 bordered minor 零と選択 minor 非零から rank と span が厳密に従い、係数は adjugate/determinant の式として指定できる。

本段階は source、所有数、正の事前 bound だけであり、55値も rank も計算しない。
実行前には独立 source/code/custody peer と ROOT GO が必要である。
full lower rank は upper2やactual HF2の等号、T2全長を決めない。
