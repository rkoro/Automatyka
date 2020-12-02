mode(0)
// Bezpieczne pobranie probki.
function [y] = pobierz_probke(x, index, default)
    try
        y = x(index);
    catch
        y = default;
    end
endfunction
// Wykrycie czasow narastania dla nastaw regulatora.
function [t10, t90] = wykrycie_narastania(syg_zad, ywyj, k, tp)
    mode(1)
    t10 = 0;
    t90 = 0;

    for i = 1:1:size(ywyj)(1)
        if ywyj(i) >= 0.9*k then t90 = i*tp; 
            break
        end
    end
    for i = 1:1:size(ywyj)(1)
        if ywyj(i) >= 0.1*k then t10 = i*tp;
            break
        end
    end
    if t10 >= t90 then t90 = 1; t10 = 1
    end
    mode(0)
endfunction
// Obliczenie nastaw zadanego regulatora lub przepisanie nastaw.
function [kp, ti, td]= dobierz_nast_reg(R, t90, t10, kp_arg, ti_arg, td_arg)
    if t10 >= t90 then 
        disp("Uklad nie osiaga 90% wartosci zadanej na tej dlugosci sygnalu")
        [kp, ti, td]= dobierz_nast_reg(R, 2, 1, kp_arg, ti_arg, td_arg)
    end
    if R == "P" then
        if kp_arg == "A" then kp = (t90) / t10; else kp = kp_arg end 
        ti = %inf;
        td = 0;
    elseif R == "PI" then
        if kp_arg == "A" then kp = 0.9 * (t90-t10) / t10 else kp = kp_arg end 
        if kp_arg == "A" then ti = t10 / 0.3; else ti = ti_arg end
        td = 0;
    elseif R == "PID" then
        if kp_arg == "A" then kp = 1.2 * (t90-t10) / t10 else kp = kp_arg; end
        if ti_arg == "A" then ti = 2 * t10 else ti = ti_arg end
        if td_arg == "A" then td = 0.5 * t10 else td = td_arg end
    else
        ti = %inf;
        td = 0;
        kp = 1;
        text = "Bledny rodzaj regulatora."
    end
    if kp == 0; then kp = 1; end
endfunction
// Tworzenie domyslnego sygnalu zadanego/wejsciowego.
function [r] = stworz_def_syg(tp, czas_pomiaru)
    j = 3.14;
    ilosc_prob = (czas_pomiaru/tp);
    for i = 1:1:ilosc_prob*0.05
        r(i) = 0;
    end

    for i = ilosc_prob*0.05:1:ilosc_prob*0.30
        r(i) = 1;
    end

    for i = ilosc_prob*0.30:1:ilosc_prob*0.7
        r(i) = cos((j*j*j*j)*(i - ilosc_prob*0.40)/(ilosc_prob*0.40)*1);
        j = j + (0.0038*tp);
    end
    j = 0;
    for i = ilosc_prob*0.7:1:ilosc_prob*0.85
        j = j + 1.2/(ilosc_prob*0.15)
        r(i) = j;
    end
    
    for i = ilosc_prob*0.85:1:ilosc_prob
        r(i) = 0;
    end
    r(ilosc_prob*0.85+ilosc_prob*0.05) = 1;
endfunction
// Wyliczenie kolejnej probki wyjsca obiektu, rownanie rekursywne.
function [yn] = obiekt (y, up, kt, a, n)
    yn = kt*up;
    for i = 1:n
        yn = yn - a(i) * y(i) 
    end
endfunction
// Wyliczenie kolejnej probki sterującej obiektem, wyjście regulatora.
function [u] = regulator (e, ep, ec, kp, ti, td, tp)
    u = kp*(e+(tp/ti)*ec+(td/tp)*(e-ep));
endfunction
// Badanie.
function [ywyj, ywyj_bez_reg, u, e] = badanie_obj(r_syg_zadany, wzmoc_dyskr, wspol_wiel_wekt, okres_probkowania, rodzaj_reg, wzmocnienie_reg, pid_ti, pid_td, czas_symulacji)
    mode(0)
    rzad_obiektu = size(wspol_wiel_wekt)(2);
    ywyj_bez_reg = [];
    ywyj = [];
    u = [];
    e = [];
    ec = 0;
    t10 = 0;
    t90 = 0;
    syg_zad = [];
    y = [];
    y=matrix(y,[0])
    J = 0;
    //Jeśli podano 0 lub "A" jako sygnal, tworzony jest sygnał zdefiniowany w funkcji, który zostanie podany na obiekt.
    if r_syg_zadany == "A" then 
        disp("Automatyczne tworzenie sygnalu zadanego ...")

        syg_zad = stworz_def_syg(okres_probkowania, czas_symulacji);
    else syg_zad = r_syg_zadany;
    end

    //Blok tworzący automatycznie nastawy regulatora P/I/D.
    if wzmocnienie_reg == "A" || pid_ti == "A" || pid_td == "A" then
        for i = 1:1:10
            r(i) = 0;
        end
        for i = 10:1:czas_symulacji/okres_probkowania
            r(i) = 1;
        end
        //Tworzenie odpowiedzi skokowej obiektu w celu wykrycia czasów narastania
        for i = 1:1:size(r)(1)
            rp = pobierz_probke(r, i-rzad_obiektu, 0);
            y = [];
            for ij = 1:1:rzad_obiektu
                y = [y, pobierz_probke(ywyj_bez_reg, i-ij, 0)];
            end 
            ywyj_bez_reg(i) = obiekt (y, rp, wzmoc_dyskr, wspol_wiel_wekt, rzad_obiektu);

    end
        // Wykrywanie czasow narastania i opcjonalne automatyczne dobranie nastaw regulatora.
        [t10, t90] =  wykrycie_narastania(r, ywyj_bez_reg, mean(ywyj_bez_reg), okres_probkowania);
        disp("Automatyczny dobor nastaw regulatora metodą Zieglera - Nicholsa ...", rodzaj_reg)

        // Jesli nastawy sa liczbami, sa przepisywane, jesli zawieraja string "A", wtedy program dobiera nastawy automatycznie.
        [kp, ti, td] = dobierz_nast_reg(rodzaj_reg, t90, t10, wzmocnienie_reg, pid_ti, pid_td)
    else
        [kp, ti, td] = dobierz_nast_reg(rodzaj_reg, 2, 1, wzmocnienie_reg, pid_ti, pid_td)
    end
        // Emulacja obiektu bez regulacji
    for i = 1:1:size(syg_zad)(1)
        rp = pobierz_probke(syg_zad, i-rzad_obiektu, 0);
        y = [];
        for ij = 1:1:rzad_obiektu
            y = [y, pobierz_probke(ywyj_bez_reg, i-ij, 0)];
        end 
        ywyj_bez_reg(i) = obiekt (y, rp, wzmoc_dyskr, wspol_wiel_wekt, rzad_obiektu);
    end
        // Emulacja obiektu z regulacja
    for i = 1:1:size(syg_zad)(1)

        up = pobierz_probke(u, i-rzad_obiektu, 0);
        y = [];
        for ij = 1:1:rzad_obiektu
            y = [y, pobierz_probke(ywyj, i-ij, 0)];
        end 

        ep = pobierz_probke(e, i-1, 0);
        ywyj(i) = obiekt (y, up, wzmoc_dyskr, wspol_wiel_wekt, rzad_obiektu);
        e(i) = syg_zad(i) - ywyj(i);
        u(i) = regulator(e(i), ep, ec ,kp, ti, td, okres_probkowania);
        ec = ec + e(i);
        J(i) = pobierz_probke(J, i-1, 0) + (e(i)^2)*okres_probkowania; //wskaznik jakosci
    end
    // Tworzenie wykresow.
    f6=scf(666);
        subplot(211);
         plot([syg_zad], 'r', 'LineWidth',1.47);
         plot([ywyj], 'b', 'LineWidth',1.47);
         xgrid;
         xlabel('Numer probki');
         ylabel('Poziom sygnalu');
         title('Sygnał wyjsciowy i wartość zadana');
         legend('Wartość zadana','sygnał wyjściowy z regulacją', 4);
        subplot(212);
         plot([u], 'r', 'LineWidth', 1.47);
         xgrid;
         xlabel('Numer probki');
         ylabel('Poziom sygnalu');
         title('Sygnał sterujący/ wyjsciowy z regulatora');
         legend('Sygnał sterujący', 4);
         
     f1=scf(111);
         subplot(211);
          plot([ywyj_bez_reg], 'r', 'LineWidth',1.47);
          plot([ywyj], 'b', 'LineWidth',1.47);
          xgrid;
          xlabel('Numer probki');
          ylabel('Poziom sygnalu');
          title('Porownanie ukladu bez regulacji i z');
          legend('Sygnał wyjściowy bez regulacji','Sygnał wyjściowy z regulacją', 4);
        subplot(212);
          plot([J], 'r', 'LineWidth', 1.47);
          plot([e], 'b', 'LineWidth', 1.47);
          xgrid;
          xlabel('Numer probki');
          ylabel('Poziom sygnalu');
          title('Wskaźnik jakości');
          legend('Wskaźnik jakości ','Uchyb', 4);
       
endfunction
//Domyslne nastawy:
disp('WARTOSCI DOMYSLNE STWORZONYCH NASTAW:')

mode(1)
syg_zadany = "A";
wzmocnienie_obiektu = 5;
T1 = 1.5;
T2 = 1.1;
okres_probkowania = 0.01;
rodzaj_reg = "PID";
pid_ti = "A";
pid_td = "A";
wzmocnienie_reg = "A";
czas_symulacji = 100;
mode(0)

r_syg_zadany = stworz_def_syg(okres_probkowania, czas_symulacji);
a = (okres_probkowania*(T1+T2) - 2*T1*T2)/(T1*T2);
b = (okres_probkowania*okres_probkowania - okres_probkowania*(T1+T2) + T1*T2)/(T1*T2);
c = 0.000009;//Współczynniki dla obiektu 5 inercyjnego, nie wywolane domyslnie
d = 0.0009;
e = 0.0000001;

kt = (wzmocnienie_obiektu*okres_probkowania*okres_probkowania)/(T1*T2);

[ywyj, ywyj_bez_reg, ster_syg, uchyb_syg] = badanie_obj(syg_zadany, kt, [a b], okres_probkowania, rodzaj_reg, wzmocnienie_reg, pid_ti, pid_td, czas_symulacji);
disp('Wywolanie:')
disp('[ywyj, ywyj_bez_reg, ster_syg, uchyb_syg] = badanie_obj(r_syg_zadany, kt, [a b], okres_probkowania, rodzaj_reg, wzmocnienie_reg, pid_ti, pid_td, czas_symulacji);')
disp('uwaga na echo')

