# DERIVATION.md — Ayrık kinetik-enerji denetimi: türetim ve konvansiyon (v0.1, 2026-09-15; Aşama 1 adım 1)

Durum: **v0.1 DOĞRULANDI (2026-09-15)** — `stage1/toyClosure.py` (2B periyodik, N = 8 düzgün / 12 ve 16 rastgele pertürbe,
rastgele alanlar): §2–§7'nin tüm özdeşlikleri (Euler/BDF2/CN kapalı formlar, midPoint/linear/upwind yüz ayrışımı, G^T
transpoze + X_f hücre özdeşliği, difüzyon non-orth dahil, tam kapanış hücre ve küresel, e_φ) rel. hata ≤ 2e-14. Adım 1 GO.
Kaynak: CLAUDE_REVIEW.md §2 iskeleti + kendi türetmelerim; OpenFOAM v2406 operatör tanımları kaynak koddan.

## 0. Notasyon
- Hücre P (hacim V_P, merkez x_P); yüz f (alan vektörü S_f owner O(f) → neighbour N(f), |S_f|, merkez x_f);
  s_{Pf} = +1 (P = O), −1 (P = N). Cyclic ve processor yüzler iç yüz gibi (komşu değerleri patchNeighbourField).
- Yüz akısı φ_f (hacimsel, owner'dan dışa pozitif) = pEqn'in ürettiği Rhie–Chow akısı. (div φ)_P V_P = Σ_f s φ_f.
- Zaman: a = U^{n+1}, b = U^n, d = U^{n−1}; Δt sabit; K_P = ½|a_P|² V_P; küresel K = Σ K_P / V_Ω.
- Linear ağırlık w_f (owner tarafı): U_f^{lin} = w_f U_O + (1−w_f) U_N; OpenFOAM: w_f = |S_f·(x_N − x_f)| / (|S_f·(x_f − x_O)| + |S_f·(x_N − x_f)|). midPoint: U_f^{mid} = ½(U_O + U_N). Düzgün hex'te w_f = ½.
- Seçilen şema: U_f^{sch} = U_f^{mid} + δ_f (δ_f: şemaya özgü sapma; limitli/gradyanlı şemalarda U'ya bağlı).

## 1. Çözücünün kurduğu ayrık momentum denklemi (pimpleFoam v2406, laminer, PISO nOuter = 1)
UEqn.H: fvm::ddt(U) + fvm::div(phi, U) + divDevReff(U) = 0, sonra solve(UEqn == −grad p); divDevReff =
−fvc::div(ν dev2(T(∇U))) − fvm::laplacian(ν, U). Hücre formunda (V_P ile çarpılmış):
  T(a)_P V_P + Σ_f s φ^n_f U_f^{sch} − Σ_f s J_f(a) − D2_P V_P + Σ_f s S_f p_f = r_P V_P          (1)
- φ_used = φ^n (UEqn dış döngü başında, pEqn'den önce kurulur; PIMPLE nOuter > 1'de son dış iterasyonun φ'si ≈ φ^{n+1}).
- Difüzyon akısı J_f = ν_f |S_f| [δ_f (U_N − U_O) + k_f·(∇U)_f]  (Gauss linear corrected; δ_f = nonOrthDeltaCoeffs,
  k_f = non-ortogonal düzeltme vektörü, düzgün hex'te 0; (∇U)_f linear interpolasyonlu Gauss gradyanı).
- D2_P = [fvc::div(ν dev2(T(∇U)))]_P (explicit; sürekli limitte ν∇(∇·u) = 0, ayrık olarak ≠ 0).
- r_P: momentum artığı — (1)'i nihai a, p^{n+1}, φ^n ile yeniden kurunca kalan. İçeriği: PISO bölme artığı (Issa 1986),
  lineer çözücü artıkları, doğrusal-olmama gecikmesi (limiter/gradyan/non-orth/dev2 terimleri çözücüde U^n veya U* ile,
  FO'da a ile değerlendirilir).
Birincil konvansiyon: (1) a_P ile nokta-çarpılır ve hacim ağırlıklı toplanır.

## 2. Zamansal terim: e_time := a·T(a) V − (K^{n+1} − K^n)/Δt
### 2.1 Euler: T(a) = (a − b)/Δt
  a·(a−b) = ½|a|² − ½|b|² + ½|a−b|²  ⇒  e_time = ½|a−b|² V/Δt ≥ 0  (kesin, disipatif).
### 2.2 backward (BDF2), sabit Δt: T(a) = (3a − 4b + d)/(2Δt)
Özdeşlik (Dahlquist G-kararlılığı; elle doğrulandı): 2a·(3a−4b+d) = |a|² + |2a−b|² − |b|² − |2b−d|² + |a−2b+d|².
E_G(a,b) := ¼(|a|² + |2a−b|²), X(a,b) := E_G(a,b) − ½|a|² = ¼(a−b)·(3a−b). Böylece
  a·T(a) = [E_G(a,b) − E_G(b,d)]/Δt + |a−2b+d|²/(4Δt)
  e_time = [X(a,b) − X(b,d)] V/Δt  +  |a−2b+d|² V/(4Δt)  =: e_time^{store} + e_time^{diss}.
- e_time^{diss} ≥ 0 kesin disipatif; e_time^{store} teleskoplar (zaman integrali sınır terimine iner), işareti belirsiz.
- Kontrol: a−b = b−d = v ⇒ e_time = ½|v|² V/Δt (yalnız depolama). ✓
- OpenFOAM değişken-adım katsayıları (ω = Δt_n/Δt_{n−1}): c = (1+2ω)/(1+ω), c0 = 1+ω, c00 = ω²/(1+ω). Sabit adımda
  3/2, 2, 1/2. Değişken adımda e_time artık olarak yine kesin hesaplanır; diss/store ayrımı ω'ya bağlı G-norm ister
  (Liao & Zhang 2021) → üretimde sabit Δt.
- İlk adım: timeIndex < 2 → deltaT0 = GREAT → c ≈ 1, c00 ≈ 0 → Euler; FO adım 1'de §2.1'i kullanır.
### 2.3 CrankNicolson(ψ): T(a) = ((1+ψ)/Δt)(a − b) − ψ·ddt0^n
ddt0^n = kayıt defterindeki `ddt0(U)` (adımın ilk fvm/fvc::ddt çağrısında güncellenir; aynı timeIndex'te tekrar
güncellenmez — evaluate() bayrağı). a·(a−b) özdeşliğiyle
  e_time = ψ[(K^{n+1}−K^n)/Δt − a·ddt0^n V] + ((1+ψ)/(2Δt))|a−b|² V.
ψ = 0 → Euler. İlk adım: coef = 1 (Euler) → §2.1.
### 2.4 Alternatif çarpan U^{n+½} = ½(a+b)
½(a+b)·(a−b)/Δt = (K^{n+1}−K^n)/Δt kesin ⇒ Euler'de e_time^{(½)} ≡ 0; zamansal hata diğer terimlere taşınır.
Yalnız küresel zaman serisi için hesaplanır (yerel yüz ayrışımı a-tabanlıdır).

## 3. Konvektif terim: W_P := a_P·Σ_f s φ_f U_f^{sch}
U_f^{sch} = U_f^{mid} + δ_f ve a_P·δ_f = ½(a_P − a_N)·δ_f + U_f^{mid}·δ_f ile, yüz başına (P = O için; N için s → −s):
  s φ_f a_P·U_f^{sch} = s φ_f ½|a_P|²  +  ½ d_f  +  s F^C_f,
  d_f := s_{Pf} φ_f (a_P − a_N)·δ_f  (P ↔ N simetrik: her iki hücrede aynı değer),
  F^C_f := φ_f [½ a_O·a_N + U_f^{mid}·δ_f]  (P ↔ N simetrik ⇒ Σ_P Σ_f s F^C_f = 0, teleskoplar).
Tanımlar (hücre başına):
  e_cont(P) := ½|a_P|² Σ_f s φ_f = ½|a_P|² (div φ)_P V_P;
  m_f := d_f[linear] = s φ_f (w_f − ½)(a_P − a_N)·(a_O − a_N) → düzgün hex'te 0;   e_mesh(P) := ½ Σ_f m_f;
  e_conv(P) := ½ Σ_f (d_f^{sch} − m_f)  (şemanın linear'dan sapması; linear için ≡ 0; midPoint için = −e_mesh).
  W_P = e_cont + e_conv + e_mesh + Σ_f s F^C_f;  küresel: Σ_P W_P = Σ_P (e_cont + e_conv + e_mesh).
Özel şemalar: upwind: δ_f = ±½(a_O − a_N) (akış yönüne göre) ⇒ d_f = ½|φ_f||a_O − a_N|² ≥ 0.
linear: δ_f = (w_f − ½)(a_O − a_N). cubic/LUST/linearUpwind/limitedLinear: δ_f gradyan/limiter içerir → d_f işaretsiz.
Hücre dağıtımı konvansiyonu: yüz terimi ½–½ (birincil); alternatif owner-tarafı (yerel alanı değiştirir, küresel toplamı
değiştirmez). Akı konvansiyonu: φ^n (birincil); e_φ(P) := a_P·Σ_f s (φ^{n+1}_f − φ^n_f) U_f^{sch}[φ^{n+1}] − a_P·Σ_f s φ^n_f U_f^{sch}[φ^n]
(genel; lineer şemalarda a_P·Σ_f s (φ^{n+1}−φ^n)_f U_f^{sch}).

## 4. Basınç terimi: a_P·(G p)_P V_P = Σ_f s S_f·a_P p_f, p_f = w_f p_O + (1−w_f) p_N
Transpoze (summation-by-parts): Σ_P a_P·(Gp)_P V_P = Σ_P p_P (G^T a)_P V_P ile
  (G^T a)_P V_P := Σ_f s_{Pf} ω_{Pf} S_f·(a_P − a_{P'}) = Σ_{f∈P} ω_{Pf} S_f·(a_O(f) − a_N(f)),   ω_{Pf} = w_f (P = O) veya 1 − w_f (P = N).
  (Uygulama notu: her iki tarafta da ω_P·S_f·(a_O − a_N) — komşu tarafında işaret ters yazılırsa özdeşlik bozulur; toy test bunu yakaladı.)
Eşdeğer: (G^T a)_P V_P = −Σ_f s S_f·Ũ_f, Ũ_f = (1−ω_{Pf}) a_P + ω_{Pf} a_{P'} (ağırlıkları yer değiştirilmiş
interpolasyon) — yani "ağırlık-transpoze ıraksama"; düzgün hex'te (w = ½) −(div_{mid} a)_P V_P. [İşaret: G^T a = −div^{dual} a.]
Hücre özdeşliği (kesin): a_P·(Gp)_P V_P = p_P (G^T a)_P V_P + Σ_f s X_f,  X_f := w_f p_O (S_f·a_N) + (1−w_f) p_N (S_f·a_O)
(P ↔ N için aynı X_f ⇒ teleskoplar).
  e_pres(P) := p_P (G^T a)_P V_P.
Fiziksel okuma: pEqn div φ = 0 sağlarken kollokasyonlu hızın ağırlık-transpoze ıraksaması sıfır değildir; fark Rhie–Chow
filtresi − ddtCorr'dur (φ = S·interp(HbyA) + ddtCorr − rAU_f ∇_f p·S; a = HbyA − rAU ∇p). ddtCorr payı ayrı ölçülmez;
"ddtCorr off" (`backward 0`) koşusuyla fark alınır. Alternatif yerel dağıtım (yüz-artığı formu): e_pres'(P) := Σ_f s p_f S_f·(a_P − U_f^{mid});
küresel toplam aynı, yerel farklı — "convention sensitivity".

## 5. Difüzyon: a_P·Σ_f s J_f
a_P·J_f = ½(a_O + a_N)·J_f + ½ s (a_O − a_N)·J_f... (P = O için s = +1). Owner: a_O·J = ½(a_O+a_N)·J + ½(a_O−a_N)·J;
neighbour: −a_N·J = −½(a_O+a_N)·J + ½(a_O−a_N)·J. Simetrik kısım F^D_f := ½(a_O + a_N)·J_f teleskoplar; her iki hücre
½(a_O − a_N)·J_f alır. J_f'nin ortogonal kısmı ν|S_f|δ_f(a_N − a_O) ile ½(a_O−a_N)·J_f = −½ g_f + ½(a_O−a_N)·ν_f|S_f| k_f·(∇a)_f,
g_f := ν_f|S_f|δ_f|a_N − a_O|² ≥ 0. Tanımlar:
  ε_ν(P) := ½ Σ_f [ g_f − (a_O − a_N)·ν_f|S_f| k_f·(∇a)_f ]   (yüz-tabanlı, implicit Laplacian ile tutarlı; poli mesh'te
  yerel ≥ 0 garanti yok — not düşülür);   a_P·Σ_f s J_f = −ε_ν(P) + Σ_f s F^D_f.
  e_diff(P) := ε_ν(P) − 2ν (S:S)_P V_P  (TANI; kapanışa girmez);  e_dev(P) := −a_P·D2_P V_P (kapanışa girer, küçük).

## 6. Artık: e_iter(P) := a_P·r_P V_P
r "çözüldüğü gibi" (φ^n, §1). Bileşenleri: (i) PISO bölme artığı; (ii) lineer çözücü artıkları (solverPerformance'tan
ayrıca raporlanır); (iii) doğrusal-olmama gecikmesi (limiter/explicit gradyan/non-orth/dev2). B0 (linear, düzgün hex,
ortogonal) için (iii) yalnız dev2 → e_iter ≈ saf bölme artığı; LUST/limitedLinear/poli'de (iii) büyür → V4 ölçer.
Adı "algebraic/iterative residual"; "dissipation" denmez.

## 7. Kapanış
(1)·a_P V_P ve §2–§6 ile, hücre başına:
  (K_P^{n+1} − K_P^n)/Δt = −e_time − e_conv − e_mesh − e_cont − ε_ν − e_dev − e_pres + e_iter − Σ_f s (F^C_f − F^D_f + X_f)
Küresel (teleskoplayan akılar sıfır):
  dK/dt = −ε_ν − e_time − e_conv − e_mesh − e_cont − e_pres − e_dev + e_iter.                                   (2)
R_closure := |LHS − RHS| / max(|LHS|, |RHS|, ε_mak); kümülatif: |K(t) − K(0) + ∫_0^t (ε_ν + Σe − e_iter) dt'| / K(0).
(2) cebirsel özdeşliktir; R_closure yüz teleskoplaması, cyclic/processor yüz işleme ve kayan-nokta indirgemesini test
eder. FO'nun çözücünün gerçekten çözdüğü denklemi gördüğü V4 (enstrümante çözücü) ile test edilir.

## 8. Konvansiyon çiftleri (aynı koşuda küresel seri; yerel alan yalnız birincil)
| Seçim | Birincil | Alternatif |
|---|---|---|
| Çarpan | U^{n+1} | U^{n+½} (§2.4) |
| Konvektif referans | midPoint (e_conv = şema − linear; e_mesh = linear − midPoint) | linear (e_mesh ≡ 0, e_conv = şema − linear) |
| Akı | φ^n | φ^{n+1} (e_φ ile) |
| Yüz→hücre dağıtımı | ½–½ | owner |
| e_pres yerel formu | p_P (G^T a)_P V_P | Σ_f s p_f S_f·(a_P − U_f^{mid}) |

## 9. İşaret ve tanım özeti (FO kodu için)
| Terim | Tanım (hücre, ×V_P) | İşaret |
|---|---|---|
| e_time | a·T(a)V − ΔK/Δt (kapalı formlar §2) | Euler ≥ 0; BDF2 diss ≥ 0, store ± |
| e_conv | ½Σ_f (d_f^{sch} − m_f) | ± (upwind: ≥ 0 toplam) |
| e_mesh | ½Σ_f m_f | ± (düzgün hex: 0) |
| e_cont | ½|a|² (div φ) V | ± (div φ ≈ p-tol) |
| e_pres | p (G^T a) V | ± |
| ε_ν | ½Σ_f [g_f − (a_O−a_N)·ν|S_f|k_f·(∇a)_f] | ≥ 0 (ortogonal); poli ± |
| e_dev | −a·D2 V | ± küçük |
| e_iter | a·r V | ± |
| e_φ | a·[C(φ^{n+1}) − C(φ^n)]a V | ± (raporlanır, (2)'ye girmez) |
| e_diff | ε_ν − 2ν(S:S)V | tanı |

## 10. Doğrulama planı (adım 1 go ölçütü)
`stage1/toyClosure.py`: 2B periyodik N×N mesh (düzgün ve rastgele pertürbe), rastgele a, b, d, p, φ, ν, ddt0; kontrol:
(i) BDF2 özdeşliği; (ii) e_time kapalı formları (Euler/BDF2/CN) ≡ tanım; (iii) konvektif yüz ayrışımı (midPoint/linear/
upwind) hücre ve küresel; (iv) basınç transpoze özdeşliği ve X_f hücre özdeşliği; (v) difüzyon hücre özdeşliği
(non-orth dahil); (vi) tam kapanış (2) ve hücre formu ≤ 1e-13 (rel.).


## Addendum (v1.1, 2026-09-18): BDF2 with a variable time step
OpenFOAM's `backward` scheme with Δt ≠ Δt₀ (r = Δt/Δt₀) uses T(a) = (c a − c₀ b + c₀₀ d)/Δt with
c = 1 + r/(1+r), c₀₀ = r²/(1+r), c₀ = c + c₀₀ (first step: Euler). The exact time term is e_time = a·T(a) V − ΔK/Δt (operator
route, verified against `fvc::ddt(U)` by `chkTime`). The constant-step split e_time = [X(a,b) − X(b,d)]V/Δt + |a − 2b + d|²V/(4Δt),
X(u,v) = ¼(u−v)·(3u−v), is an identity only for r = 1; for r ≠ 1 the function object keeps the G-norm change as the storage part
and reports the remainder e_time − storage as the dissipative part (Grigorieff: sign-definiteness of variable-step BDF2 requires a
step-dependent G-norm, which would not telescope). Test: 16³ TGV, `adjustTimeStep yes, maxCo 0.5`, Δt varying 0.12 → 0.22 over
7 steps: R_closure ≤ 3e-15 at every step (v1.0: 5e-2 at the steps where Δt changed).
