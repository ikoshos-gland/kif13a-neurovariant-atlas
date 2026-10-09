# KIF13A K1122E: hesaplamalar ve MD durumu

4 Ekim 2026. Gerçek RNA ve ThermoMPNN tahminleri TRUBA hesaplama düğümlerinde çalıştırıldı, kaynak hash'leri doğrulanarak yerel depoya alındı. Bağımsız deneysel panel analizi tamamlandı. Uzun MD koşuları gönderildi; henüz bunlardan tamamlanmış bilimsel sonuç yok.

## 1. Doğru transkriptte RNA etkisi

Referans MANE Select **ENST00000259711.11 / NM_022113.6**, GRCh38, negatif zincir. **chr6:17785639 T>C = c.3364A>G = AAA→GAA = K1122E** eşlemesi doğrulandı. Değişim exon 28'in üçüncü bazında; bu konum splice etkisini kendiliğinden kanıtlamaz. 1121. kodon iki exon arasında bölünüyor; düz bir genomik offset varsayılmadı.

SpliceAI'nin beş ve Pangolin v2'nin 12 resmi ağırlığı, sabitlenmiş kaynak sürümlerinden indirildi. Modeller doğru genomik bağlamın transkripsiyon yönünde çalıştırıldı; çıktılar genomik koordinatlara geri çevrildi. ±50 ve ±500 baz sonuçları, ham ve exon sınırlarıyla maskelenmiş skorlar ve tüm konumsal profiller kaydedildi. Modellerin farklı maskeleme kuralları birbirine karıştırılmadı.

**K1122E:**

- SpliceAI maksimum delta: **0,005326**, iki aralıkta da aynı.
- Pangolin maksimum ham değişim: **0,003538**. Maskeli maksimum: ±50 bazda **0,002249**, ±500 bazda **0,003520**.
- Pangolin beyin kanalında maksimum ham değişim yaklaşık **0,000992**. Bu kanal belirli bir beyin hücre tipinin modeli değildir.

Sekiz girdi kullanıldı: WT→WT, K1122E, K1122R, eş anlamlı K1122K, popülasyonda gözlenen K1104E/E1128K ve exon 28'in kanonik acceptor/donor bozma kontrolleri. WT→WT sıfır verdi; kanonik bozmalarda SpliceAI yaklaşık **0,99**, Pangolin kayıp skorları **0,86** düzeyinde. Bunlar yayımlanmış hasta RNA ölçümleri olmayan sentetik sistem kontrolleridir, klinik duyarlılık/özgüllük hesabı değildir. Eş anlamlı K1122K için SpliceAI **0,2493**, Pangolin **0,0673** çıktı; eş anlamlılık RNA bakımından negatif kontrol sayılmadı.

**Yorum:** İki model K1122E için güçlü bir splice değişimini desteklemiyor. Hasta RNA'sında küçük, hücreye özgü veya modellerin kaçırdığı etkiler dışlanmış değil. Kaynaklar: [SpliceAI](https://github.com/illumina/SpliceAI), [Pangolin](https://github.com/tkzeng/Pangolin).

## 2. Bölgesel karşılaştırma

1100-1138 arasında **39 × 19 = 741** aminoasit değişimi ve **351** olası kodlayıcı SNV üretildi. Çok baz değişikliği gerektiren aminoasitler de protein tablosunda var; her aminoasit değişiminin tek bir DNA değişimine karşılık geldiği varsayılmadı.

AlphaMissense tüm 741 değişime eklendi. K1122E skoru **0,8044**, bölgede yüksekten düşüğe **249/741**; bu bölge için doğrulanmış klinik kesim değeri yok. Mevcut VEP önbelleğinden 55 değişime ESM1b eklendi; **K1122E için ESM1b yok**, eksik değer sıfırla doldurulmadı.

203 Compara ortolog hizalamasında insan MANE kaynağı kontrol edildi. 1122. pozisyonda standart aminoasit taşıyan 191 dizinin **186'sı K**, ikisi R, biri E. K oranı **%97,38**; kalan 12 giriş gap/belirsiz aminoasit. Alternatife özgü sözde sayımlı log-frekans oranları üretildi. Bu diziler filogenetik olarak bağımsız gözlemler değildir; oran klinik olasılık veya kalibre edilmiş etki skoru değildir.

## 3. Protein stabilitesi ve etkileşim ayrımı

Resmi ThermoMPNN ağırlıklarıyla, AlphaFold v6 monomerinin **1092-1146** ve **1084-1154** fragmanlarında aynı 741 değişim hesaplandı: **1.482** tahmin. Girdi dizisi/kalıntı eşlemesi ve WT→WT sıfır kontrolü geçti. Aşağıdakiler modelin kcal/mol ölçeğindeki stabilite değişimi tahminleridir; pozitif işaret destabilizasyon yönündedir.

| Değişim | AlphaMissense | ΔΔG ana fragman | ΔΔG geniş fragman | gnomAD gözlenen AC |
|---|---:|---:|---:|---:|
| K1122E | 0,8044 | +1,645 | +1,671 | gözlenmedi |
| K1122R | 0,1926 | +0,137 | +0,132 | gözlenmedi |
| K1104E | 0,3609 | +0,284 | +0,343 | 1 |
| E1128K | 0,2333 | −0,069 | −0,095 | 27 |

K1122E ana fragmanda **82/741**, aynı pozisyonda **2/19** yüksek ΔΔG sırasındadır. Yükü koruyan K1122R'nin tahmini çok daha küçük. Dolayısıyla bulgu yalnız “her yük değişimi eşittir” biçiminde açıklanamaz; konum önem taşıyor. İki fragmanın K1122E farkı **0,0256 kcal/mol**. Bu benzerlik, fragman varsayımının veya modelin hedefte doğru olduğunu kanıtlamaz.

Model, izole tahmini fragmandaki stabiliteyi değerlendirir; tam proteinin katlanma enerjisi veya partner bağlanma serbest enerjisi ölçülmedi. Yapısal AF modelinde K1122 NZ ile E1130 karboksilat oksijenleri yaklaşık **2,65 Å** yakınlıkta; bu yalnız statik, monomer içi tuz köprüsü hipotezidir. MD'nin önceden tanımlı temas tabloları bu çiftin korunmasını da ölçüyor. K1122 çevresinde doğrulanmış bir cargo/dimer geometrisi olmadığı için bağlanma kaybı veya allosteri kanıtı çıkarılmadı. Yüksek etki skoru ve düşük stabilite skoru tek başına allosteri demek değildir. [ThermoMPNN yöntem/kod](https://github.com/Kuhlman-Lab/ThermoMPNN), [resmi skor işareti dönüşümü](https://github.com/Kuhlman-Lab/ThermoMPNN/blob/370f76ec62bd929f7425e311d8df04a0d094990f/datasets.py).

## 4. Gerçek MD hazırlığı ve üretim

Ana model kullanıcı talebiyle **ff19SB+OPC** olarak güncellendi: OpenMM **8.3.1**, PDBFixer **1.12**, MDTraj **1.11**; 310 K, 1 bar, 0,15 M tuz, 2 fs adım. OpenMM 8.3.1 seçilen Amber19 dosyalarını içeriyor. ACE/NME nötr kaplar ff19SB ile beş değişim × iki sınır için on gerçek yerel denetimden ve on CUDA denetiminden sonlu enerjiyle geçti. Hidrojen yükleyicisinin NME tanımını ezme sorunu korunmuş düzeltmeyle gideriliyor. OPC yükleri, dört noktalı su yapısı, kütlesiz sanal noktaları ve ff19SB CMAP omurga terimleri sistem kurulurken denetleniyor. Başlangıç su koordinatları TIP4P-Ew kutusundan; gerçek su/iyon parametreleri `amber19/opc.xml` dosyasından geliyor. Bu koordinat şablonu fiziksel su modelini değiştirmiyor. Parametre dosyalarının hash'leri kaydediliyor. [OpenMM resmi kuvvet alanı belgeleri](https://docs.openmm.org/latest/userguide/application/02_running_sims.html).

Güncel **v3** protokolünde, serbestçe dönen fragmanın periyodik kopyalarla yapay temasını azaltmak için bütün eksenlerde en uzun başlangıç atom boyutu + 2,4 nm olan kübik kutu korunuyor. Önceki ff14SB/TIP3P v2 protokolü ve pilotu tarihsel kayıt olarak ayrı tutuluyor. AF güveni core ortalaması **83,47**, K1122 **80,31**; geniş fragmanın uçlarında **49,41** kadar düşük pLDDT var, dolayısıyla o sınır duyarlılık kontrolüdür. Sistem izole kuyruk fragmanı, su ve iyonlardan oluşuyor; mikrotübül, motor bölgesi, cargo partneri veya dimer içermiyor.

Seçimin gerekçesi aminoasit türüne özgü omurga enerjileri ve heliks eğilimleridir; K1122E için yerel heliks/yük değişimi sorusuyla uyumu yöntemsel bir çıkarımdır. ff19SB yöntem çalışması OPC birlikteliğini öneriyor. 2026'da yüklü helikal peptitleri inceleyen bağımsız çalışma hem ff14SB/TIP3P hem ff19SB/OPC'nin protonlanmaya bağlı değişimleri yakalayabildiğini, heliks miktarı ve hidrojen bağı kararlılığında farklılıklar bulunduğunu bildiriyor. Bu, KIF13A'da ff19SB'nin üstünlüğünü kanıtlamaz. [ff19SB yöntem çalışması](https://pmc.ncbi.nlm.nih.gov/articles/PMC13071887/), [2026 helikal peptit karşılaştırması](https://pubmed.ncbi.nlm.nih.gov/42267453/). Başlangıç pH 7 protonlanması sabit tutuluyor; constant-pH MD yapılmıyor.

WT, K1122E, K1122R, K1104E ve E1128K için iki fragmanda üç bağımsız tohum: **30 koşu**. Ana fragman **100 ns**, geniş fragman **50 ns**; toplam **2,25 µs** üretim hedefi. Üretim öncesi 1,2 ns dengeleme; ilk 10 ns analiz dışı. Aynı AF omurgası başlangıcından bağımsız hız/iyon/hidrojen örneklemesi kullanılıyor; farklı katlanma havzalarının örneklendiği varsayılmıyor.

Birincil ölçümler core alfa torsiyon oranı ve i→i+4 hidrojen bağı oranı. İkincil ölçümler yerel heliks oranı, ağır atom temasları, tuz köprüleri, CA RMSD. İki sınırın tekrar karşılaştırmaları ortak **10-50 ns** penceresini kullanıyor; bütün sürelerin ortalamaları/temas tabloları ayrıca betimleyici çıktı olarak korunuyor. Otokorelasyon, bloklar ve tekrarlar arası farklar kaydediliyor. Güven aralıkları **tekrar ortalamalarından** kuruluyor; kareler bağımsız örnek sayılmıyor. Periyodik öz temas veya fragmanlar arasında tutarsız etki, sağlam mekanizma yorumunu engelliyor. Üç tekrar keşifsel kesinlik sağlar; toplam potansiyeller mutant-WT arasında katlanma/bağlanma ΔΔG'si olarak çıkarılmıyor. Eski ff14SB/GBn2 statik enerji tanılaması v3 ana modelinde kapalı; ThermoMPNN ayrı stabilite tahmini olarak korunuyor.

TRUBA ff19SB/OPC pilot dizisi **6444622**: ana WT, ana K1122E ve geniş WT, her biri 1,2 ns tam dengeleme + 0,1 ns üretim. Üçü de CUDA üzerinde **COMPLETED 0:0**, üç pilot denetimi geçti. Parçacık sayıları sırasıyla **58.009, 58.036 ve 94.676**; son raporlanan hızlar **208, 204 ve 128 ns/gün**. Pilot denetimi sonlu enerji/kayıtlar, 290-330 K sıcaklık, pozitif hacim, gerçek OPC parametreleri ve ff19SB omurga terimlerini gerektiriyor. Bu pilotlar bilimsel çıkarımdan dışlanıyor. Üretim dizisi **6444625** için pilot bağımlılığı sağlandı; analiz/arşiv **6444626**, bütün üretimin başarısına bağlı. İş başına süre yeni model için güvenlik payıyla **36 saat**. Eski **6443760/6443761** işleri üretim başlamadan iptal edildi. Önceki başarılı **6443669** pilotu ff14SB/TIP3P kullandığından onun 277 ns/gün hızı OPC için kullanılamaz.

Kuyruk bilgisi tarih damgalı `data/interim/k1122e_md/TRUBA_jobs_20261004.json` dosyasında. Son kontrolde üç GPU pilotu tamamlandı, üretim **Priority**, analiz **Dependency** nedeniyle bekliyor. Tamamlanmış pilotların küçük kayıt arşivi `data/raw/completed_opc_pilots_6444622.tar.gz` olarak alındı. Eski başlangıç tahmini yeni dizi için geçerli değil. Çalışma: `/arf/scratch/avural/kif13a_k1122e_20261004`; v3 kodu `scripts/k1122e_md_v3`, çıktıları ayrı `_v3` dizinlerinde. Biten MD sonuçları tek arşivle `/arf/home/avural/kif13a_k1122e_20261004_release` altında korunacak. Üretim sonuçları henüz tamamlanmadı.

## 5. MD'nin ek değerinin bağımsız sınanması

avGFP, Bgl3, GB1, Pab1 ve Ube4b için **7.511** yayımlanmış tek değişim işlev ölçümü, yazarların Rosetta özellikleri ve gerçek yayımlanmış WT-MD RMSF dosyaları kullanıldı. Bütün kalıntılar ve ayrıca bütün proteinler dışarıda bırakılarak test yapıldı. Ölçekleme/eğitim yalnız eğitimde, hiperparametre seçimi yalnız doğrulama kümesinde; tüm kalıntıyla eşleştirilmiş bootstrap kullanıldı. Kaynak değişimler sıfır tabanlı, RMSF'nin “residue” sütunu atom kimliğiydi; ikisi doğrulanarak düzeltildi.

Enerji içeren temel modelde, kalıntı testindeki Spearman Bgl3 için **0,427→0,462**, Pab1 için **0,240→0,294**; avGFP ve GB1 iyileşmedi. Gerçek düşük maliyetli aminoasit özellikleri temelinde Bgl3 **0,243→0,319**, Pab1 **0,130→0,245**. Sonuncusu, ilk sonuçlar görüldükten sonra tanımlanmış ve ayrıca beyan edilmiş keşifsel duyarlılık analizidir. Daha geniş %99,5 aralıklar, ilk analizde Ube4b'nin yeni-protein olumlu sonucunu doğrulamıyor. Ek katkı bağlama ve temelin gücüne bağlı.

Bu sınama WT esnekliğinin ek bilgisini değerlendiriyor; K1122E mutant MD'sinin doğrulanması, hedef kuyruğa aktarım veya modern dil modellerinden üstünlük testi değildir. [Yayımlanmış yöntem](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1013728), [yazar verileri](https://github.com/SBarethiya/VEP). Ayrıntılar `data/interim/k1122e_benchmark/` altındaki iki raporda.

## Çıkarılabilen sonuç

Tamamlanan hesaplar K1122E için RNA etkisinden çok **aday bir yerel protein stabilitesi/elektrostatik etkiyi** destekliyor. Bu, model temelli hipotezdir. Bitmemiş MD'den sonuç çıkarılmadı; hedefin deneysel işlev/abundance/binding ölçümleri yok. HudsonAlpha kaydındaki VUS, kamuya açık bilgilerle klinik neden ilan edilmedi; başka varyant veya nedenler dışlanmış değil.

`scripts/k1122e_rna_region/`, `scripts/k1122e_energy/`, `scripts/k1122e_md/`, `scripts/k1122e_benchmark/` ve `scripts/truba/k1122e/` yürütülebilir betikleri içeriyor. Girdiler/çıktılar kaynak ve SHA256 yan dosyalarını taşıyor; birleşik 741 satır ve kontrol tablosu `data/interim/k1122e_integrated/`, şekil `figures/k1122e/` altında.

## Yeniden üretim ve doğrulama

Yerel veri hazırlığı: `python -m scripts.k1122e_rna_region.prepare`. Modellerin gerçek çıkarımı `scripts/truba/k1122e/03_inference.slurm` ile SLURM üzerinde yapıldı. Sonuçların birleştirilmesi: `python -m scripts.k1122e_energy.integrate`. Bağımsız panel: `python -m scripts.k1122e_benchmark.fetch`, ardından aynı paketin `run`, `report`, `low_cost_sensitivity` ve `low_cost_report` modülleri. Panel girdileri `data/raw/k1122e_benchmark/source` altında, hash manifesti sonuç klasöründe.

V3 denetimi: `python -m unittest discover -s tests -p 'test_k1122e_*.py' -v` ile 24 test, **22 başarılı**, OpenMM gerektiren iki test atlandı. Aynı MD test dosyası izole OpenMM ortamında ayrıca **10/10 başarılı**: gerçek AlphaFold ff19SB kapları, OPC sanal noktaları ve sonlu kısa dinamik, yanlış TIP4P-Ew parametrelerinin reddi ve ortak zaman penceresi dahil. Ayrıca on tam varyant/sınır kap denetimi yerelde ve CUDA'da geçti. Güncel hash/kayıt `data/interim/k1122e_md_v3/validation_20261004.json`; önceki v2 denetimi `data/interim/k1122e_integrated/validation_20261004.json` dosyasında tarihsel olarak korunuyor. Kalıntı ayrımı, eğitim dışı sızıntı, negatif zincir splice koordinatları ve pilotların bilimsel analizden dışlanması da denetlendi.

Yayında TRUBA kaynağı belirtilmeli: “Bu araştırmada yer alan tüm/kısmi nümerik hesaplamalar TÜBİTAK ULAKBİM, Yüksek Başarım ve Grid Hesaplama Merkezi'nde (TRUBA kaynaklarında) gerçekleştirilmiştir.”
