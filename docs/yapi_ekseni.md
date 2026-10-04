# Yapı ekseni: KIF13A motor alanındaki adayların üç boyutlu konumu

Tarih: 16 Eyl 2026. Üretici: `scripts/structure_axis.py`. Çıktılar:
`data/interim/structure/KIF13A_motor_residue_structure_20260916.tsv` (kalıntı başına),
`data/interim/structure/KIF13A_candidates_structure_20260916.tsv` (aday başına),
`data/interim/structure/KIF13A_structure_axis_summary_20260916.json` (çalışma özeti),
`figures/KIF13A_structure_axis.svg`.

Amaç: `data/processed/KIF13A_candidates_20260915.tsv` düz bir sıralı listeydi. Bir varyantın motor alanın
neresinde durduğunu söylemiyordu. Bu eksen her adaya bir mekanizma etiketi ekler: ATP cebi mi, mikrotübül
arayüzü mü, gömülü çekirdek mi, yoksa yüzey mi.

## 1. Neden ödünç bir yapı

KIF13A motor alanının deneysel yapısı yok (`data/raw/pdb/` içinde Q9H1H9 için kapsam dosyası bile yok).
Paraloğu KIF1A'nın ise mikrotübüle bağlı kriyo-EM yapıları var. Bu yüzden geometri KIF1A'dan ödünç alınır:
KIF13A AlphaFold modeli, KIF1A zincirinin üstüne bindirilir ve mesafeler o çerçevede ölçülür.

Kullanılan yapı: **8UTN**, `KIF1A[1-393] AMP-PNP bound two-heads-bound state in complex with a microtubule
(class T23L1)`, 3,1 Å kriyo-EM. Birincil atıf, indirilen mmCIF dosyasının kendi `_citation` bloğundan alındı
(uydurulmadı) ve Europe PMC üzerinden canlı doğrulandı: Benoit MPMH, Rao L, Asenjo AB, Gennerich A, Sosa H.
Nat Commun 2024;15:5530, doi `10.1038/s41467-024-48720-4`, PMID `38956021`.

Zincirler dosyanın `_struct_ref` / `_entity` bloklarından okunur, varsayılmaz:

| Zincir | Varlık | UniProt |
|---|---|---|
| K, N | Kinesin-like protein KIF1A | Q12756 (kalıntı 1-393) |
| A, E, S | Tubulin alpha-1B chain | Q2XVP4 (domuz) |
| B, I | Tubulin beta-2B chain | A0A287AGU7 (domuz) |

Heteroatomlar: K zincirinde ANP (AMP-PNP) + Mg, N zincirinde ANP; tübülinlerde GTP/GDP/taksol.
**Not:** 8UTX KIF1A değildir, kullanılmadı. KIF1A 1-393 numaralandırması MANE (Q12756-3) ile birebir aynıdır
(`CLAUDE.md`), yani paralog haritası doğrudan kullanılabilir.

## 2. Yöntem

1. **Bindirme.** `data/interim/paralog_maps/KIF1A_to_KIF13A_pairwise_MANE.tsv` haritasındaki, KIF1A tarafı
   1-365 aralığına düşen ve her iki yapıda da CA atomu bulunan kalıntı çiftleri alınır. En küçük kareler
   bindirmesi yapılır, sonra sapması 3 Å'dan büyük çiftler atılıp bindirme yinelenir (ıraksak ilmekler
   çekirdek hizalamasını bozmasın diye).
2. **Baş çerçeveleri.** İki başlı durumda K ve N başlarının ikisi de gerçekten mikrotübüle bağlıdır ve
   tübüle biraz farklı kalıntılarla dokunurlar. Bu yüzden nükleotid mesafesi yalnız birincil çerçeveden
   (K, çünkü Mg + ANP orada ve uyumu en iyi), tübül mesafesi ise **iki çerçevenin en küçüğünden** alınır.
   Her çerçevenin kendi değeri de tabloda ayrı sütun olarak durur (`dist_tubulin_chainK_A`,
   `dist_tubulin_chainN_A`).
3. **Mesafeler.** Her KIF13A kalıntısı için, ağır atomlar üzerinden nükleotide (ANP + Mg) ve herhangi bir
   tübülin atomuna en kısa uzaklık.
4. **Gömülülük.** Shrake-Rupley (Bio.PDB.SASA) **yalıtılmış motor alan** üzerinde çalıştırılır (kalıntı
   1-380). Tüm 1805 aa'lık model kullanılsaydı, düşük güvenli kuyruk motorun üstüne serilip sahte gömülülük
   yaratırdı. Bağıl SASA için Tien MZ ve ark., PLoS One 2013;8:e80635 azami ASA değerleri kullanıldı
   (doi `10.1371/journal.pone.0080635`, canlı doğrulandı). Yedek ölçü olarak 10 Å
   içindeki CB komşu sayısı da yazılır.
5. **Güven.** Kalıntı pLDDT değeri AlphaFold güven dosyasından
   (`AF-Q9H1H9-F1-confidence_v6.json`) okunur, B-faktöründen türetilmez.

### Bindirme kalitesi

| Çerçeve | Çift sayısı | RMSD (ham) | Budanmış çift | RMSD (budanmış) |
|---|---|---|---|---|
| K zinciri (birincil) | 354 | 2,55 Å | 311 | **1,03 Å** |
| N zinciri | 356 | 8,85 Å | 228 | 1,70 Å |

354 çiftin 221'i özdeş kalıntı. Budamada atılan 43 çift ağırlıklı olarak bilinen ıraksak ilmeklerde
toplanıyor (KIF13A numarasıyla 42-46, 236-242, 257-269, 296-300, 361-363); yani budama rastgele değil,
gerçekten farklı olan bölgeleri ayıklıyor. N zincirinin ham RMSD'sinin yüksek olması beklenen bir şey:
öndeki baş gergin bir konformasyonda ve boyun bağlayıcısı yerleşmemiş durumda, o yüzden daha az çift hayatta
kalıyor. Her iki çerçeve de "birkaç Ångström" ölçütünü geçtiği için ikisi de kullanıldı; `--max-rmsd`
(varsayılan 3 Å) geçilseydi betik hesap yapmadan durur.

### Eşikler

| Sınıf | Kural |
|---|---|
| ATP cebi | nükleotide (ANP + Mg) en kısa ağır-atom uzaklığı ≤ 5,0 Å |
| mikrotübül arayüzü | tübüline en kısa ağır-atom uzaklığı ≤ 5,0 Å (yukarıdaki değilse) |
| gömülü çekirdek | bağıl SASA < 0,15 (yukarıdakilerin ikisi de değilse) |
| yüzey | geri kalanlar |
| boyun bağlayıcı (güvenilmez) | pozisyon > 352, aşağıdaki kısıtlara bakınız |

5 Å, temas için alışılmış bir kesme değeri: hidrojen bağı ve tuz köprüsü mesafelerini (2,5-3,5 Å) ve van der
Waals temasını (3,5-4,5 Å) içine alır, ikinci kabuk su aracılı temasları dışarıda bırakır. 0,15 bağıl SASA
gömülülük için standart eşiktir. Sıralama önemlidir: tek bir kalıntı (KIF13A 257) hem ATP hem tübül eşiğini
geçiyor, ATP olarak sınıflandırıldı. Ham sayılar tabloda olduğu için isteyen başka eşikle yeniden
sınıflandırabilir.

## 3. Sağlama kontrolleri

Aynı ölçüm hem KIF13A'nın bindirilmiş kalıntısı için hem de **yapının içindeki gerçek KIF1A kalıntısı** için
yapıldı. Böylece bir tuhaflığın haritalamadan mı, bindirmeden mi, yoksa gerçek biyolojiden mi geldiği
ayırt edilebiliyor.

| KIF1A | KIF13A | Gerçek KIF1A ölçümü | Bindirilmiş KIF13A ölçümü | Sonuç |
|---|---|---|---|---|
| T99 | T104 | ATP 3,0 Å (K başı) | ATP 2,1 Å | **ATP cebi**, beklendiği gibi |
| R254 | R258 | MT 2,9 Å (α-tübülin E420) | MT 0,6 Å | **mikrotübül arayüzü**, ama çakışma var (aşağıya bakınız) |
| P305 | P303 | MT 5,3 Å (K), 6,7 Å (N) | MT 4,9 Å | **mikrotübül arayüzü**, sınırda |
| R350 | R348 | MT 7,4 Å (K başı), **2,4 Å (N başı, α-tübülin E415)** | MT 1,4 Å | **mikrotübül arayüzü** |
| P364 | P362 | MT 14,8 Å, ATP 39,2 Å | MT 12,6 Å | boyun bağlayıcı, arayüzde değil, doğru |

Dördü de beklenen yere düştü. En öğretici olanı R350: literatürdeki kinezin-tübülin tuz köprüsü iddiasını
**yapının kendisinde doğruladım**, ama yalnız N başında. K başında aynı arginin tübülden 7,4 Å uzakta.
Yani tuz köprüsü baş durumuna bağlı ve tek bir baş üzerinden ölçüm yapılsaydı R348 kaçırılırdı. İki
çerçeveli ölçümün gerekçesi budur. Karşı taraf α-tübülin E415'tir (arginin-glutamat, kimyasal olarak
anlamlı); en yakın atom sıralamasında E415'i α-tübülin R402 izler, o bir tuz köprüsü değildir.

P362 (KIF1A P364) ATP cebinden de arayüzden de uzak çıkıyor. Bu doğru: boyun bağlayıcı prolini mekanizmayı
temas üzerinden değil, boyun bağlayıcının yerleşmesi üzerinden etkiler. Yapı ekseni bu mekanizmayı
ölçemez, o yüzden 352'den sonrası ayrı sınıfa konuldu.

### Ölçülen ek doğrulamalar

- 8UTN'deki gerçek KIF1A başı tübüle dokunuyor: K başında 30, N başında 36 kalıntı 5 Å içinde;
  1-352 aralığında iki başın birleşimi 38 kalıntı.
- Bindirilmiş KIF13A'da motor alan içinde 50 kalıntı arayüzde çıkıyor. Gerçek ayak izinden (38) yaklaşık
  %30 fazla. Fark, en iyimser yorumla bile şişme payıdır; kaynağı aşağıdaki çakışmalar.
- 17 KIF13A kalıntısı tübüline 2 Å'dan yakın düşüyor (258-260, 269-270, 273, 298-301, 344, 348, 376-380).
  Bunlar fiziksel değil: katı-cisim yerleştirmede yan zincir rotamerleri yeniden optimize edilmiyor.
  Mesafeyi "temas var/yok" olarak okumak doğru, "0,6 Å" sayısını bir bağ uzunluğu gibi okumak yanlış.

## 4. Sonuçlar

Motor alan (5-352, 348 kalıntı): ATP cebi 25, mikrotübül arayüzü 50, gömülü çekirdek 124, yüzey 149.
Motor alanın pLDDT ortancası 90,0; 53 kalıntı 70'in altında.

Adayların dağılımı (parantez içi: gnomAD'de güvenilir gözlenenler):

| Katman | n | ATP cebi | MT arayüzü | Gömülü çekirdek | Yüzey | Boyun / motor dışı |
|---|---|---|---|---|---|---|
| A | 52 (26) | 23 (7) | 18 (11) | 9 (6) | 1 (1) | 1 (1) |
| B | 73 (43) | 5 (4) | 8 (6) | 37 (19) | 20 (11) | 3 (3) |
| C | 170 (110) | 1 (1) | 5 (5) | 4 (2) | 12 (4) | 148 (98) |
| D | 33 (23) | 3 (2) | 2 (0) | 21 (15) | 7 (6) | 0 |

Okunuşu:

- **Katman A gerçekten mekanik.** 52 adayın 41'i (%79) ya ATP cebinde ya da mikrotübül arayüzünde. Yalnız
  bir tanesi düpedüz yüzeyde: M214R, ve onun pLDDT'si 42, yani zaten güvenilmez bir bölge. Bu, A katmanının
  (AlphaMissense + paralog uyumu) bağımsız bir kanıt ekseninde de tutarlı olduğunu gösteriyor.
- **Katman B farklı bir mekanizma sınıfı.** B'de baskın sınıf gömülü çekirdek (%51) ve yüzey (%27).
  Yani B, katlanma kararlılığını bozan ya da işlevi belirsiz olan değişimlerden oluşuyor; A ise doğrudan
  katalitik ve arayüz kalıntılarından. İkisi aynı listenin iki ucu değil, iki ayrı hipotez sınıfı.
- **Katman C zaten büyük ölçüde motor alan dışında** (170'in 148'i), beklendiği gibi.
- **Katman D** (paralog kanıtı var, AlphaMissense eşiğin altında) çoğunlukla gömülü çekirdek. Modelin eşiği
  geçmemesi ile kalıntının gömülü olması birbiriyle çelişmiyor, ama D'nin A'ya terfi etmesi için yapısal
  kanıt yeterli değil.

### Katman A adayları

`evet/hayır` sütunu gnomAD'de güvenilir gözlem (AC0 filtresi yok) anlamına gelir. Mesafeler Å.

| Varyant | Yapısal sınıf | ATP | Tübül (en yakın) | Bağıl SASA | pLDDT | Güvenilir | KIF1A'da patojenik |
|---|---|---|---|---|---|---|---|
| R11G | ATP cebi | 3,1 | 15,5 | 0,04 | 97 | hayır | R11L, R11Q, R11W |
| R11L | ATP cebi | 3,1 | 15,5 | 0,04 | 97 | hayır | R11L, R11Q, R11W |
| R11P | ATP cebi | 3,1 | 15,5 | 0,04 | 97 | hayır | R11L, R11Q, R11W |
| R11Q | ATP cebi | 3,1 | 15,5 | 0,04 | 97 | hayır | R11L, R11Q, R11W |
| R11W | ATP cebi | 3,1 | 15,5 | 0,04 | 97 | **evet** | R11L, R11Q, R11W |
| R13G | ATP cebi | 3,1 | 15,1 | 0,15 | 94 | hayır | R13L, R13C, R13S, R13H |
| R13L | ATP cebi | 3,1 | 15,1 | 0,15 | 94 | **evet** | R13L, R13C, R13S, R13H |
| R13P | ATP cebi | 3,1 | 15,1 | 0,15 | 94 | hayır | R13L, R13C, R13S, R13H |
| R13Q | ATP cebi | 3,1 | 15,1 | 0,15 | 94 | hayır | R13L, R13C, R13S, R13H |
| R13W | ATP cebi | 3,1 | 15,1 | 0,15 | 94 | hayır | R13L, R13C, R13S, R13H |
| P14H | ATP cebi | 3,0 | 19,2 | 0,37 | 91 | hayır | P14A, P14T |
| P14L | ATP cebi | 3,0 | 19,2 | 0,37 | 91 | hayır | P14A, P14T |
| P14S | ATP cebi | 3,0 | 19,2 | 0,37 | 91 | hayır | P14A, P14T |
| P14T | ATP cebi | 3,0 | 19,2 | 0,37 | 91 | hayır | P14A, P14T |
| E19A | gömülü çekirdek | 8,0 | 15,7 | 0,05 | 86 | **evet** | E19K |
| E19K | gömülü çekirdek | 8,0 | 15,7 | 0,05 | 86 | hayır | E19K |
| S63P | ATP cebi | 3,2 | 21,8 | 0,06 | 90 | hayır | S58W, S58L |
| S63Y | ATP cebi | 3,2 | 21,8 | 0,06 | 90 | hayır | S58W, S58L |
| A73P | gömülü çekirdek | 5,6 | 24,2 | 0,03 | 90 | **evet** | A68E |
| Q75H | gömülü çekirdek | 5,8 | 23,4 | 0,04 | 93 | **evet** | Q70K |
| T104K | ATP cebi | 2,1 | 9,1 | 0,07 | 83 | hayır | T99K, T99M |
| G107V | ATP cebi | 2,8 | 14,7 | 0,16 | 90 | hayır | G102S, G102D |
| R173G | mikrotübül arayüzü | 24,8 | 2,4 (β E410) | 0,41 | 87 | **evet** | R169W, R169T, R169K |
| R173L | mikrotübül arayüzü | 24,8 | 2,4 (β E410) | 0,41 | 87 | hayır | R169W, R169T, R169K |
| R173P | mikrotübül arayüzü | 24,8 | 2,4 (β E410) | 0,41 | 87 | hayır | R169W, R169T, R169K |
| H175L | mikrotübül arayüzü | 29,2 | 3,0 (β M406) | 0,16 | 76 | hayır | H171P |
| H175R | mikrotübül arayüzü | 29,2 | 3,0 (β M406) | 0,16 | 76 | **evet** | H171P |
| M214R | yüzey | 5,4 | 8,6 | 0,28 | **42** | **evet** | M210R |
| S219N | ATP cebi | 2,1 | 11,5 | 0,10 | 62 | **evet** | S215N, S215R |
| R220C | ATP cebi | 4,3 | 7,3 | 0,11 | 66 | **evet** | R216H, R216C |
| R220H | ATP cebi | 4,3 | 7,3 | 0,11 | 66 | **evet** | R216H, R216C |
| H222R | gömülü çekirdek | 7,8 | 11,2 | 0,02 | 87 | **evet** | H218R |
| A254V | ATP cebi | 4,3 | 9,7 | 0,00 | 86 | **evet** | A250P, A250S, A250D |
| E257K | ATP cebi | 4,9 | 2,6 | 0,09 | 75 | **evet** | E253D, E253K |
| R258S | mikrotübül arayüzü | 10,9 | 0,6 (α G416) | 0,62 | 67 | **evet** | R254G, R254P, R254Q, R254W |
| G272S | mikrotübül arayüzü | 11,8 | 4,3 (α V409) | 0,02 | 64 | **evet** | G268W, G268R, G268E |
| S278L | gömülü çekirdek | 12,3 | 6,5 | 0,02 | 85 | **evet** | S274L |
| S278W | gömülü çekirdek | 12,3 | 6,5 | 0,02 | 85 | hayır | S274L |
| L279I | gömülü çekirdek | 8,2 | 6,0 | 0,05 | 86 | hayır | L275Q, L275P |
| G283V | mikrotübül arayüzü | 15,7 | 3,6 (α E415) | 0,07 | 90 | **evet** | G279R |
| G283W | mikrotübül arayüzü | 15,7 | 3,6 (α E415) | 0,07 | 90 | hayır | G279R |
| P303R | mikrotübül arayüzü | 21,5 | 4,9 (β F260) | 0,08 | 90 | **evet** | P305L |
| P303S | mikrotübül arayüzü | 21,5 | 4,9 (β F260) | 0,08 | 90 | hayır | P305L |
| P303T | mikrotübül arayüzü | 21,5 | 4,9 (β F260) | 0,08 | 90 | hayır | P305L |
| R305L | mikrotübül arayüzü | 22,6 | 2,4 (β N414) | 0,56 | 92 | hayır | R307P, R307G, R307Q |
| R305Q | mikrotübül arayüzü | 22,6 | 2,4 (β N414) | 0,56 | 92 | **evet** | R307P, R307G, R307Q |
| D306H | mikrotübül arayüzü | 20,3 | 2,1 (β P261) | 0,40 | 90 | **evet** | D308G |
| D306N | mikrotübül arayüzü | 20,3 | 2,1 (β P261) | 0,40 | 90 | **evet** | D308G |
| T310I | gömülü çekirdek | 17,0 | 8,8 | 0,01 | 94 | **evet** | T312I, T312P |
| R348G | mikrotübül arayüzü | 19,4 | 1,4 (α R402) | 0,37 | 93 | **evet** | R350W, R350G |
| R348Q | mikrotübül arayüzü | 19,4 | 1,4 (α R402) | 0,37 | 93 | **evet** | R350W, R350G |
| L370P | boyun bağlayıcı (güvenilmez) | 35,1 | 4,8 | 0,34 | 79 | **evet** | L372V |

Dikkat çeken üç nokta:

1. **R11, R13, P14 kümesi.** KIF13A'nın N ucundaki bu üç kalıntı nükleotide 3 Å uzaklıkta çıkıyor. Bu bir
   bindirme yanılgısı değil: 8UTN'deki gerçek KIF1A R11, R13, P14 kalıntıları da AMP-PNP'ye 3,5-3,8 Å
   uzaklıkta. Yani kinezin motorunun N uç şeridi gerçekten nükleotid cebinin kenarında duruyor. Toplam 14
   katman A adayı bu üç kalıntıda toplanıyor ama yalnız ikisi (R11W, R13L) gnomAD'de güvenilir biçimde
   gözlenmiş, geri kalanı AC0 filtreli.
2. **R348 (KIF1A R350).** 2026 Nisan'ında yayımlanan R350G ve R350W yapıları (9YA5, 9YA7, 9YAB, 9YAI)
   bu kalıntının önemini gösteriyor. KIF13A'daki karşılığı R348G ve R348Q hem güvenilir gözlenmiş hem de
   arayüzde. Bu, listedeki en doğrudan test edilebilir hipotez.
3. **Düşük pLDDT kümesi 214-222 ve 258-272.** S219N, R220C/H, R258S, G272S gibi adaylar pLDDT 62-67
   bandında. Yapısal etiketleri doğru yönü gösteriyor olabilir ama tek başlarına kanıt sayılmamalı;
   `structure_evidence_ok` sütunu bunları zaten `False` işaretliyor.

## 5. Kısıtlar

- **Bu bir KIF13A yapısı değil.** KIF1A çerçevesine oturtulmuş bir KIF13A modelidir. Motor çekirdek
  1,03 Å'ya oturuyor ama ıraksak ilmeklerde (42-46, 236-242, 257-269, 296-300) hata birkaç Ångström.
  O ilmeklerdeki mesafeleri sınıf sınırına yakınsa (4-6 Å) kararsız kabul etmek gerekir.
- **Yan zincirler yeniden optimize edilmedi.** Katı-cisim yerleştirme 17 kalıntıda 2 Å'nın altında çakışma
  üretiyor ve arayüz sayısını gerçek KIF1A ayak izine göre %30 şişiriyor. "Temas var" kararı sağlam,
  mesafenin tam değeri değil.
- **Tek bir nükleotid durumu.** 8UTN AMP-PNP bağlı; ATP cebi mesafeleri bu duruma özgü. Apo (8UTS) ve ADP
  (8UTR) yapılarında cep açılır. Betik `--pdb-id` ile o yapılarla da çalıştırılabilir, henüz çalıştırılmadı.
- **Boyun bağlayıcı ölçülemiyor.** 353-380 arası AlphaFold'da tek bir konformasyonda; kriyo-EM'de yerleşmiş
  (K) ve yerleşmemiş (N) iki farklı yol izliyor. Bu yüzden o bölge ayrı sınıfa kondu. L370P (katman A)
  bu nedenle "mikrotübül arayüzünde" diye rapor edilmiyor, hesaplanan 4,8 Å tabloda duruyor ama iddiaya
  çevrilmiyor.
- **Tübülin izotipi domuzdan.** 8UTN'de α-1B (Q2XVP4) ve β-2B (A0A287AGU7) domuz tübülini var. Arayüz
  kalıntıları türler arasında yüksek derecede korunmuştur ama en yakın tübülin kalıntısının numarası
  insan izotipinde birebir aynı olmayabilir.
- **Gömülülük yalıtılmış motor alandan.** Gerçek KIF13A'da motor alan, boyun ve CC1 bölgesiyle komşudur;
  gerçek gömülülük burada hesaplanandan biraz yüksek olabilir.
- **Bu eksen bir patojenite skoru değil.** Hiçbir sınıf tek başına "zararlı" demek değildir; A katmanındaki
  adayların çoğunun işlevsel bir bölgede olması, o eksenin bağımsız bir doğrulaması olarak okunmalıdır.

## 6. Yeniden üretme

```
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\structure_axis.py
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\structure_axis.py --pdb-id 8UTS --kif1a-chain K
mcp\atlas_server\.venv\Scripts\python.exe -u scripts\structure_axis.py --atp-cutoff 4.5 --mt-cutoff 4.5
```

8UTN.cif bir kez `data/raw/pdb/` içine iner ve provenans kaydı `data/raw/PROVENANCE.jsonl` dosyasına
yazılır; sonraki çalışmalar yeniden indirmez. Bütün çıktı tabloları `atlas_mcp.core.save_tsv` ile
yazıldığı için yanlarında `.prov.json` dosyaları vardır.
