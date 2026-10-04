# Kalibrasyon ölçümleri: KIF1A setleri üzerinde model performansı

Tarih: 15 Eyl 2026. Setler: `data/interim/clinvar_curation/KIF1A_calibration_sets_20260915.tsv`
(kürasyon: `docs/KIF1A_ClinVar_kurasyonu.md`). Skorlar: tek bir Ensembl VEP çalışmasıyla, bütün varyantlar
için aynı kaynaktan (`KIF1A_calibration_scores_20260915.tsv`). Metrikler:
`KIF1A_calibration_metrics_20260915.json`.

Yöntem notları:
- Güven aralıkları **kalıntı-gruplu bootstrap** (2.000 tekrar): aynı kalıntının farklı değişimleri birlikte
  yeniden örneklenir, yani tekrarlayan kalıntılar (ör. R254) yapay güven yaratmaz.
- `sens@90%spec`: popülasyon kontrollerinde %90 özgüllük veren eşikte duyarlılık. Eşik, o modelin bu gen
  için kullanılabilir kesme değeridir.
- ESM1b skoru ters yönlüdür (daha negatif = daha zararlı), hesapta yönü çevrildi.

## 1. Ana ölçüm: motor domain içi

Pozitif: ClinVar çekirdek P/LP missense (motor, n=110). Negatif: gnomAD v4.1 motor missense, AC ≥ 3 (n=61).

| Model | AUROC | %95 GA | sens@90%spec | Eşik |
|---|---|---|---|---|
| MetaRNN (etiketle eğitilmiş) | 0,949 | 0,913-0,977 | 0,85 | 0,928 |
| AlphaMissense | **0,925** | 0,878-0,964 | 0,79 | **0,968** |
| REVEL (etiketle eğitilmiş) | 0,920 | 0,868-0,959 | 0,71 | 0,855 |
| ESM1b | 0,851 | 0,774-0,917 | 0,45 | -14,4 |
| PrimateAI | 0,801 | 0,720-0,875 | 0,37 | 0,883 |
| CADD v1.7 | 0,781 | 0,696-0,858 | 0,55 | 31 |
| Referans: KIF13A'da aynı kalıntı | 0,744 | | | |
| Evo 2 40B (16 Eyl 2026 eklendi) | 0,742 | 0,643-0,839 | 0,32 | -22,7 |

AC ≥ 2 negatifleriyle (n=93) tablo pratikte aynı: AlphaMissense 0,915, CADD 0,795, ESM1b 0,850.

## 2. Alt kümeler

| Karşılaştırma | n+ | AlphaMissense | CADD | ESM1b |
|---|---|---|---|---|
| Fonksiyonel literatürdeki pozitifler | 14 | 0,966 | 0,842 | 0,806 |
| PP3'e dayalı, ayrı tutulan pozitifler | 17 | 0,961 | 0,814 | 0,871 |
| Çekirdek pozitifler (referans) | 110 | 0,925 | 0,781 | 0,851 |

PP3'e dayanan (sınıflandırmasında in silico kanıt kullanılmış) varyantlarda AlphaMissense'in daha yüksek
çıkması beklenen döngüsellik işaretidir: bu varyantlar zaten tahmin araçları "zararlı" dediği için
sınıflandırılmış. Ana metrikte bu 17 varyant yok, doğru olan da bu.

## 3. Domain karışıklığının etkisi

Tüm protein üzerinde, ClinVar B/LB negatifleriyle:

| Model / referans | AUROC |
|---|---|
| Referans: "motor domaininde mi" | **0,995** |
| REVEL | 0,993 |
| MetaRNN | 0,986 |
| AlphaMissense | 0,959 |
| ESM1b | 0,945 |
| CADD | 0,886 |

Hiçbir model konumun tek başına verdiği ayrımı geçemiyor. Bu tablo yalnızca "neden böyle ölçmüyoruz"un
kanıtı olarak duruyor.

## 4. Konsensüs

Etiketle eğitilmemiş dört modelin (AlphaMissense, ESM1b, PrimateAI, CADD) sıra ortalaması:
**AUROC 0,911** (%95 GA 0,859-0,954), sens@90%spec 0,68. Yani eşit ağırlıklı konsensüs tek başına
AlphaMissense'i (0,925) geçmiyor; zayıf üyeler sinyali seyreltiyor. Ağırlıklandırma denenecekse
kalıntı-gruplu çapraz doğrulama içinde öğrenilmeli.

Model uyumu (Spearman, ana karşılaştırma): AlphaMissense ile MetaRNN 0,80, AlphaMissense ile REVEL 0,71,
AlphaMissense ile CADD 0,49. CADD diğerlerinden en ayrık model, ama aynı zamanda en zayıfı.

## 4b. Evo 2 40B (16 Eyl 2026)

Aynı 171 varyant, genin okunma yönünde 8.192 bazlık pencere, skor varyant bazından itibaren toplanan
log-olabilirlik farkı (`scripts/evo2_calibration.py`, çıktı `data/interim/evo2/KIF1A_calibration_evo2_w8192.tsv`).

| Ölçüm | AUROC | %95 GA |
|---|---|---|
| Evo 2, site + downstream | 0,742 | 0,643-0,839 |
| Evo 2, yalnız varyant bazı terimi | 0,761 | 0,674-0,843 |
| AlphaMissense (aynı varyantlar) | 0,925 | 0,878-0,964 |
| AlphaMissense + Evo 2 sıra ortalaması | 0,887 | 0,833-0,932 |

Sonuç: Evo 2 bu missense setinde AlphaMissense'in belirgin şekilde altında ve karıştırınca performansı
düşürüyor. AlphaMissense'in kararsız kaldığı bölgede (skor 0,50-0,99) Evo 2'nin ayrım gücü 0,611, yani
tamamlayıcı da değil. İki model arasındaki uyum Spearman 0,48.

Evo 2'ye göre en bozucu pozitifler R307P, R307G, P305L, T312P; en bozucu bulduğu negatifler N361I, D23H,
Q121P. Model prolin sokan değişimlere ve kodon bağlamına tepki veriyor gibi görünüyor, ama patojeniteyi
ayırmada yeterli değil.

Not: bu sonuç Evo 2'yi projeden çıkarmaz. Modelin iddiası kodlamayan varyantlarda; missense burada onun
güçlü olduğu alan değil. Asıl sınavı `splice_positive` seti (26 varyant) ile AF eşleşmiş intronik
negatifler üzerinde olacak; orada AlphaMissense hiç oynayamıyor.

Maliyet ve hız: 171 varyant = 342 NIM çağrısı. İlk çalıştırmada varyant başına Ensembl'den ayrı pencere
çekilince 18 dakika/varyant oldu; bölge tek seferde indirilip yerelde kesilince ve 4 çağrı paralel
yapılınca 3,3 saniye/varyanta düştü.

## 5. Zor örnekler

- **En zor pozitifler** (konsensüs en düşük): T35A, L249M, I27F, Y89F, N272S, A68E.
- **En zor negatifler** (konsensüs en yüksek, popülasyonda görülmesine rağmen): S5L, R18W, M147T, R153H,
  A206V, A9V. Bunların AF değerleri 2e-6 ile 5e-6 arasında, yani çok nadir; bir kısmı gerçekten zararlı
  olabilir (etiket gürültüsü). AC ≥ 5 seti bu riski azaltır ama n düşer.

## 6. Projeye etkileri

1. **Gen-özgü eşik şart.** AlphaMissense'in genel "likely pathogenic" sınırı 0,564'tür; bu gende %90
   özgüllük için eşik 0,968 çıkıyor. KIF13A'ya aktarımda genel sınır kullanılırsa yanlış pozitif seli olur.
2. **Konsensüse REVEL ve MetaRNN girmiyor** (etiketle eğitilmişler); ana metrikte yalnız referans olarak durur.
3. **Beklenen tavan:** bu genin motor domaininde ulaşılabilir ayrım AUROC 0,92 civarı. Plandaki 0,85 hedefi
   motor içi ölçümle de karşılanıyor, ama artık hedefin "kolay" versiyonu (0,99) ile karıştırılmayacak.
4. **Paralog korunumu tek başına 0,744.** KIF13A'ya aktarımda korunum sinyali modelden bağımsız katkı verir,
   ama modeli değerlendirirken referans çizgisi olarak raporlanmalı.

## 7. Düzeltme

`docs/KIF1A_ClinVar_kurasyonu.md` ilk sürümünde motor içi CADD AUROC 0,914 yazıyordu. O hesapta pozitiflerin
CADD skoru VEP'ten, negatiflerin ki gnomAD'in kendi anotasyonundan geliyordu; iki farklı kaynak, iki farklı
sürüm. Bütün varyantlar tek VEP çalışmasından skorlanınca CADD 0,781'e düşüyor. Ders: bir karşılaştırmanın iki
sınıfı asla farklı skor kaynaklarından beslenmemeli.

## 8. Sırada

- Evo 2 40B ile aynı setin ölçülmesi (varyant başına 2 NIM çağrısı, 171 varyant için yaklaşık 342 çağrı).
- Splice pozitifleri (26) için AF eşleşmiş intronik negatif setin kurulması, sonra SpliceAI ile Evo 2 kıyası.
- Ağırlıklı konsensüsün kalıntı-gruplu çapraz doğrulamayla öğrenilmesi.
