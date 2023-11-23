PROGRAM cas2

    IMPLICIT NONE

    !DECLARATION DE VARIABLES - GAUSS!
    DOUBLE PRECISION, ALLOCATABLE :: ag(:,:), x(:)! Augumented matrix and solution vector
    
    !DECLARATION DES VARIABLES!
    INTEGER :: nbr_noeud_x,nbr_noeud_z
    DOUBLE PRECISION :: epaisseur, longueur, teneur_eau_ini, teneur_eau_sat, masse_volumique,alpha
    DOUBLE PRECISION :: vecteur_gravite, masse_molaire_vap, cte_gaz_parfait, vitesse_air_ambiant,weq
    DOUBLE PRECISION :: temperature_air_ambiant, hr_air_ambiant,Dzz,pas_temps ,beta,temperature_humide

    DOUBLE PRECISION, EXTERNAL :: calcul_flux
    DOUBLE PRECISION, EXTERNAL :: calcul_th
  
    DOUBLE PRECISION, DIMENSION(:,:) , ALLOCATABLE :: matrice_calcul_KM, KM_l, KM_u
    DOUBLE PRECISION, DIMENSION(:) , ALLOCATABLE :: matrice_b,second_membre,LU_y
    DOUBLE PRECISION:: pas_x,pas_z,x_k,x_l,x_m,z_k,z_l,z_m,det_JE
    DOUBLE PRECISION :: flux_k,flux_l,flux_m

    INTEGER :: i,j,w,nombre_triangle,paritee,balayage_triangle,balayage_sommet,choix_methode
    INTEGER :: k,l,m,entree,test_choix_thumide,n
    INTEGER :: choix_affichage_final,choix_affichage,temps_ecriture,temps_experience,temps_x
    DOUBLE PRECISION :: a,b,c
    
    INTEGER, DIMENSION(:,:) , ALLOCATABLE :: geometrie
    DOUBLE PRECISION, DIMENSION(:,:) , ALLOCATABLE :: coord
    INTEGER, DIMENSION(:,:) , ALLOCATABLE :: connexe

    DOUBLE PRECISION, DIMENSION(:,:) , ALLOCATABLE :: matrice_masse
    DOUBLE PRECISION, DIMENSION(:,:) , ALLOCATABLE :: matrice_rigidite

    DOUBLE PRECISION, DIMENSION(:) , ALLOCATABLE :: matrice_U_N
    DOUBLE PRECISION, DIMENSION(:) , ALLOCATABLE :: matrice_U_N1
    
    !-----------------------------------------------------Conditions Initiales!-----------------------------------------------------!
    CALL lecture1(epaisseur,alpha, teneur_eau_ini, teneur_eau_sat, masse_volumique)
    CALL lecture2(vecteur_gravite,masse_molaire_vap, cte_gaz_parfait, vitesse_air_ambiant)
    CALL lecture3(temperature_air_ambiant, hr_air_ambiant,Dzz,pas_temps,temps_experience,beta,&
    &nbr_noeud_x,nbr_noeud_z,temperature_humide,temps_ecriture)
    longueur=alpha*epaisseur
    temperature_air_ambiant=temperature_air_ambiant + 273.15 !Conversion en K
    a = (-2.86d-5)*(temperature_air_ambiant**2)-(1.07d-2)*temperature_air_ambiant+10.24
    b = (-5.41d-4)*temperature_air_ambiant+1.01
    c = (4.97d-6)*temperature_air_ambiant**2-(2.67d-3)*temperature_air_ambiant+0.35
    temperature_air_ambiant=temperature_air_ambiant - 273.15 !Retour en °C
    weq=log( (-log(hr_air_ambiant)/A) - (C/-A) )/(100*log(B))

    ALLOCATE(geometrie(nbr_noeud_z,nbr_noeud_x))
    ALLOCATE(coord(nbr_noeud_x*nbr_noeud_z,3))
    ALLOCATE(connexe(2*(nbr_noeud_x-1)*(nbr_noeud_z-1),7))
    
    ALLOCATE(matrice_masse(nbr_noeud_z*nbr_noeud_x,nbr_noeud_z*nbr_noeud_x))
    ALLOCATE(matrice_rigidite(nbr_noeud_z*nbr_noeud_x,nbr_noeud_z*nbr_noeud_x))
    ALLOCATE(matrice_b(nbr_noeud_z*nbr_noeud_x))

    ALLOCATE(matrice_calcul_KM(nbr_noeud_z*nbr_noeud_x,nbr_noeud_z*nbr_noeud_x))
    ALLOCATE(KM_l(nbr_noeud_z*nbr_noeud_x,nbr_noeud_z*nbr_noeud_x))
    ALLOCATE(KM_u(nbr_noeud_z*nbr_noeud_x,nbr_noeud_z*nbr_noeud_x))
    ALLOCATE(LU_y(nbr_noeud_z*nbr_noeud_x))

    ALLOCATE(matrice_U_N(nbr_noeud_z*nbr_noeud_x))
    ALLOCATE(matrice_U_N1(nbr_noeud_z*nbr_noeud_x))
    ALLOCATE(second_membre(nbr_noeud_z*nbr_noeud_x))

    !ALLOCATE GAUSS! (AX=B)
    ALLOCATE(ag(nbr_noeud_z*nbr_noeud_x, nbr_noeud_z*nbr_noeud_x+1)) ! +1 car matrice augmenté => AX=B => MATRICE(A RESOUDRE) X = MATRCICE(SECOND MEMBRE)
    ALLOCATE(x(nbr_noeud_z*nbr_noeud_x)) !Vecteur solution X
    
    CALL execute_command_line('clear')

    PRINT*,'-------------------------------------------------------------------------------------'
    PRINT*,"                            Calculs Scientifiques II : "
    PRINT*,"   SIMULATION DU TRANSPORT D'HUMIDITE EN MILIEU POREUX LORS D'UN SECHAGE CONVECTIF"
    PRINT*,'-------------------------------------------------------------------------------------'
    PRINT*,"                     Codé par Ibrahim FSEIL et Erwan BENHENOU"
    PRINT*,"                               ENSGTI - PROMO 2024 "
    PRINT*,""
    PRINT*,"CLIQUER ENTRER"
    READ(*,*)

    CALL execute_command_line('clear')
    PRINT*,'      ---------------------------------------------------------------------'
    PRINT*,"                           Détermination de T_{humide} "
    PRINT*,'      ---------------------------------------------------------------------'
    PRINT*,"      TAPER 1 >>>>>>>>>>>>>>>>> LECTURE DANS LE FICHIER : 'donnee.txt' "
    PRINT*,"      TAPER 2 >>>>>>>>>>>>>>>>> CALCUL PAR DICHOTOMIE"
89  READ(*,*)test_choix_thumide

    IF (test_choix_thumide == 2) THEN
        a=0.0 !Borne inférieur
        b=temperature_air_ambiant !Borne superieur 
        temperature_humide= (b+a)/2 !Initilisation de th au mileu de A-B
        CALL dichotomie_Th(calcul_th,temperature_humide,a,b,hr_air_ambiant,temperature_air_ambiant)
    ELSEIF (test_choix_thumide == 1) THEN
        ! Température lu dans le fichier
    ELSE
        PRINT*,""
        PRINT*, " ERREUR, CHOISISSER 1 OU 2 !!"
        GO TO 89
        
    ENDIF

    !CALL execute_command_line('clear')
   ! PRINT*,'      ---------------------------------------------------------------------'
   ! PRINT*,"                        Méthode de resolution du système ?"
   ! PRINT*,'      ---------------------------------------------------------------------'
   ! PRINT*,"      TAPER 1 >>>>>>>>>>>>>>>>> GAUSS "
   ! PRINT*,"      TAPER 2 >>>>>>>>>>>>>>>>> DÉCOMPOSITION LU "
   ! PRINT*,"      TAPER 3 >>>>>>>>>>>>>>>>> CHOLESKY "
   ! READ(*,*)choix_methode

    choix_methode=1 !!LU et CHOLESKY ne fonctionne pas....

    !CALL execute_command_line('clear')

    PRINT*,'      ---------------------------------------------------------------------'
    PRINT*,"             Affichage du système toutes les",temps_ecriture,"secondes ?"
    PRINT*,'      ---------------------------------------------------------------------'
    PRINT*,"        OUI ========>>>>>> 1"
    PRINT*,"        NON ========>>>>>> 2"
121 READ(*,*)choix_affichage

    PRINT*,""
    PRINT*, "                              CHOIX VALIDÉ !"
    IF (choix_methode == 1) THEN
        PRINT*,""
        PRINT*, " T_{humide} =",temperature_humide,"°C et résolution par Méthode de GAUSS "
    ELSEIF (choix_methode ==2)THEN
        PRINT*,""
        PRINT*, " T_{humide} =",temperature_humide,"°C et Résolution par Méthode LU "
    ELSEIF (choix_methode ==3)THEN
        PRINT*,""
        PRINT*, " T_{humide} =",temperature_humide,"°C et Résolution par CHOLESKY "
    ELSE
        PRINT*,""
        PRINT*, " ERREUR, CHOISISSER 1,2 ou 3 !!"
        GO TO 121
    ENDIF
    PRINT*,""
    PRINT*,"CLIQUER ENTRER"
    READ(*,*)

    !------------------------------------------------------------Création du système (Mailles)------------------------------------------------------------
    !Création de coord (Tableau des coordonées des Mailles)
    !pas de la longueur!
    pas_x = longueur/(nbr_noeud_x-1)
    pas_z = epaisseur/(nbr_noeud_z-1)

    k=0
    DO j =1,nbr_noeud_x
        DO i =1,nbr_noeud_z
            k=k+1
            geometrie(i,j)=k
            coord(k,1)=k
            !Coordonnée selon x (x_i)
            coord(k,2)=(j-1)*pas_x
            !Coordonnée selon z (z_i)
            coord(k,3)=(i-1)*pas_z
        ENDDO
    ENDDO

    !------------------------------------------------------------Maillage (Création des triangles)----------------------------------------------------------
    !Création des triangles dans le premier sens
    k=0
    nombre_triangle=0
    DO i = 1,nbr_noeud_z-1 ! -1 car trangle ne peut pas depasser le tableau
        DO j = 1,nbr_noeud_x-1
            k=k+1
            nombre_triangle=nombre_triangle+1
            paritee=(i+j)/2
            connexe(k,1)=k
            IF ((i==1) .and. (j==1)) THEN !Triangle en HAUT à GAUCHE
                connexe(k,2)=geometrie(i,j)
                connexe(k,3)=geometrie(i+1,j+1)
                connexe(k,4)=geometrie(i,j+1)
            ELSEIF ((j==nbr_noeud_x-1) .and. (i==1)) THEN !Triangle en BAS à DROITE
                connexe(k,2)=geometrie(i+1,j)
                connexe(k,3)=geometrie(i,j+1)
                connexe(k,4)=geometrie(i,j)
            ELSEIF ((i==nbr_noeud_z-1) .and. (j==1)) THEN !Triangle en HAUT à GAUCHE
                connexe(k,2)=geometrie(i+1,j)
                connexe(k,3)=geometrie(i,j+1)
                connexe(k,4)=geometrie(i,j)
            ELSEIF ((i==nbr_noeud_z-1) .and. (j==nbr_noeud_x-1)) THEN !Triangle en HAUT à DROITE
                connexe(k,2)=geometrie(i,j)
                connexe(k,3)=geometrie(i+1,j+1)
                connexe(k,4)=geometrie(i+1,j)
            ELSEIF ( (i+j) == 2*paritee ) THEN !Somme ligne colone paire 
                connexe(k,2)=geometrie(i,j)
                connexe(k,3)=geometrie(i+1,j+1)
                connexe(k,4)=geometrie(i,j+1)
            ELSE    !si i+j n'est pas paire et pas un Triangle des coins => somme impaire
                connexe(k,2)=geometrie(i+1,j)
                connexe(k,3)=geometrie(i,j+1)
                connexe(k,4)=geometrie(i,j)
            ENDIF
        ENDDO
    ENDDO
    !Passage du maillage des triangles dans l'autre sens
    !k garde la valeur d'avant pour continuer la numérotation des triangles 
    DO i = 1,nbr_noeud_z-1
        DO j = 1,nbr_noeud_x-1
            k=k+1
            nombre_triangle=nombre_triangle+1
            paritee=(i+j)/2 !Permet de vérifier si i+j pair
            connexe(k,1)=nombre_triangle
            IF ((i==1) .and. (j==1)) THEN !Triangle en HAUT à GAUCHE
                connexe(k,2)=geometrie(i,j)
                connexe(k,3)=geometrie(i+1,j)
                connexe(k,4)=geometrie(i+1,j+1)
            ELSEIF ((j==nbr_noeud_x-1) .and. (i==1)) THEN !Triangle en bas à GAUCHE
                connexe(k,2)=geometrie(i+1,j)
                connexe(k,3)=geometrie(i+1,j+1)
                connexe(k,4)=geometrie(i,j+1)
            ELSEIF ((i==nbr_noeud_z-1) .and. (j==1)) THEN !Triangle en HAUT à droite
                connexe(k,2)=geometrie(i+1,j)
                connexe(k,3)=geometrie(i+1,j+1)
                connexe(k,4)=geometrie(i,j+1)
            ELSEIF ((i==nbr_noeud_z-1) .and. (j==nbr_noeud_x-1)) THEN !Triangle en bas à droite
                connexe(k,2)=geometrie(i,j)
                connexe(k,3)=geometrie(i,j+1)
                connexe(k,4)=geometrie(i+1,j+1)
            ELSEIF ( (i+j) == 2*paritee ) THEN !somme ligne colone pair
                connexe(k,2)=geometrie(i,j)
                connexe(k,3)=geometrie(i+1,j)
                connexe(k,4)=geometrie(i+1,j+1)
            ELSE    !si ce n'est pas paire et pas un Triangle des coins => somme impair 
                connexe(k,2)=geometrie(i+1,j)
                connexe(k,3)=geometrie(i+1,j+1)
                connexe(k,4)=geometrie(i,j+1)
            ENDIF
        ENDDO
    ENDDO
    
    !Création dans connexe de la présence d'un sommet sur un bord!
    DO balayage_triangle=1,nombre_triangle !Boucle sur tout les triangles
        DO balayage_sommet=2,4 !BOUCLE SUR CHAQUE SOMMET DES TRIANGLES -- car connexe(1,x)= numéro du Triangle, on part de la colonne 2
            IF (balayage_sommet==2) THEN !si c'est le sommet k
                IF (coord(connexe(balayage_triangle,2),2) == 0)THEN !x du Triangle == 0.
                    connexe(balayage_triangle,5)=1
                ELSEIF (coord(connexe(balayage_triangle,2),2) == longueur )THEN !x du Triangle == longueur
                    connexe(balayage_triangle,5)=1
                ELSEIF (coord(connexe(balayage_triangle,2),3) == 0)THEN !z du Triangle == 0.
                    connexe(balayage_triangle,5)=1
                ELSEIF (coord(connexe(balayage_triangle,2),3) == epaisseur  )THEN !z du Triangle == epaisseur
                    connexe(balayage_triangle,5)=1
                ELSE
                    connexe(balayage_triangle,5)=0
                ENDIF
            ELSEIF (balayage_sommet==3) THEN !si c'est le sommet l
                IF (coord(connexe(balayage_triangle,3),2) == 0)THEN !x du Triangle == 0.
                    connexe(balayage_triangle,6)=1
                ELSEIF (coord(connexe(balayage_triangle,3),2) == longueur )THEN !x du Triangle == longueur
                    connexe(balayage_triangle,6)=1
                ELSEIF (coord(connexe(balayage_triangle,3),3) == 0)THEN !z du Triangle == 0.
                    connexe(balayage_triangle,6)=1
                ELSEIF (coord(connexe(balayage_triangle,3),3) == epaisseur )THEN !z du Triangle == epaisseur
                    connexe(balayage_triangle,6)=1
                ELSE
                    connexe(balayage_triangle,6)=0
                ENDIF
            ELSEIF (balayage_sommet==4) THEN !si c'est le sommet m
                IF (coord(connexe(balayage_triangle,4),2) == 0)THEN 
                    connexe(balayage_triangle,7)=1
                ELSEIF (coord(connexe(balayage_triangle,4),2) == longueur )THEN 
                    connexe(balayage_triangle,7)=1
                ELSEIF (coord(connexe(balayage_triangle,4),3) == 0)THEN 
                    connexe(balayage_triangle,7)=1
                ELSEIF (coord(connexe(balayage_triangle,4),3) == epaisseur )THEN 
                    connexe(balayage_triangle,7)=1
                ELSE
                    connexe(balayage_triangle,7)=0
                ENDIF
            ENDIF
        ENDDO
    ENDDO

    !------------------------------------------------------------Création des matrices de resolution----------------------------------------------------------
    Matrice_Masse = 0
    Matrice_Rigidite = 0
    matrice_b = 0
    
    DO  balayage_triangle = 1,nombre_triangle !Avec nombre_triangle = nombre de Triangle créé
        k = connexe(balayage_triangle,2)
        l = connexe(balayage_triangle,3)
        m = connexe(balayage_triangle,4)
        x_k = coord(k,2)
        z_k = coord(k,3)
        x_l = coord(l,2)
        z_l = coord(l,3)
        x_m = coord(m,2)
        z_m = coord(m,3)
        det_JE = abs((x_l-x_k)*(z_m-z_k)-(x_m-x_k)*(z_l-z_k))
        
        Matrice_Masse(k,k)= Matrice_Masse(k,k) + (2*det_JE/24)
        Matrice_Masse(l,l)= Matrice_Masse(l,l) + (2*det_JE/24)
        Matrice_Masse(m,m)= Matrice_Masse(m,m) + (2*det_JE/24)

        Matrice_Masse(k,l)= Matrice_Masse(k,l) + (1*det_JE/24)
        Matrice_Masse(k,m)= Matrice_Masse(k,m) + (1*det_JE/24)
        Matrice_Masse(m,k)= Matrice_Masse(m,k) + (1*det_JE/24)
        Matrice_Masse(m,l)= Matrice_Masse(m,l) + (1*det_JE/24)
        Matrice_Masse(l,k)= Matrice_Masse(l,k) + (1*det_JE/24)
        Matrice_Masse(l,m)= Matrice_Masse(l,m) + (1*det_JE/24)
    
        Matrice_Rigidite(k,k)= Matrice_Rigidite(k,k) + ((pas_temps * Dzz )/(2*det_JE))&
        &*(beta*(z_l-z_m)**2+(x_m-x_l)**2)
        Matrice_Rigidite(l,l)= Matrice_Rigidite(l,l) + ((pas_temps * Dzz)/(2*det_JE))&
        &*(beta*(z_m-z_k)**2+(x_k-x_m)**2)
        Matrice_Rigidite(m,m)= Matrice_Rigidite(m,m) + ((pas_temps * Dzz)/(2*det_JE))&
        &*(beta*(z_k-z_l)**2+(x_l-x_k)**2)

        Matrice_Rigidite(k,l)= Matrice_Rigidite(k,l) + ((pas_temps *Dzz)/(2*det_JE))&
        &*(beta*(z_l-z_m)*(z_m-z_k)+(x_m-x_l)*(x_k-x_m))
        Matrice_Rigidite(k,m)= Matrice_Rigidite(k,m) + ((pas_temps *Dzz)/(2*det_JE))&
        &*(beta*(z_l-z_m)*(z_k-z_l)+(x_m-x_l)*(x_l-x_k))
        Matrice_Rigidite(m,k)= Matrice_Rigidite(m,k) + ((pas_temps *Dzz)/(2*det_JE))&
        &*(beta*(z_l-z_m)*(z_k-z_l)+(x_m-x_l)*(x_l-x_k))
        Matrice_Rigidite(m,l)= Matrice_Rigidite(m,l) + ((pas_temps *Dzz)/(2*det_JE))&
        &*(beta*(z_m-z_k)*(z_k-z_l)+(x_k-x_m)*(x_l-x_k))
        Matrice_Rigidite(l,k)= Matrice_Rigidite(l,k) + ((pas_temps *Dzz)/(2*det_JE))&
        &*(beta*(z_l-z_m)*(z_m-z_k)+(x_m-x_l)*(x_k-x_m))
        Matrice_Rigidite(l,m)= Matrice_Rigidite(l,m) + ((pas_temps *Dzz)/(2*det_JE))&
        &*(beta*(z_m-z_k)*(z_k-z_l)+(x_k-x_m)*(x_l-x_k))
    ENDDO



    !------------------------------------------------------------RESOLUTION DE (M+K)Un+1 = M Un +B----------------------------------------------------------
    !Initalisation!
    matrice_U_N=0
    DO i=1,nbr_noeud_x*nbr_noeud_z
        matrice_U_N(i)=teneur_eau_ini !à t=0, w(t=O)=w_ini en tout point
    ENDDO

    matrice_calcul_KM=0
    DO i=1,nbr_noeud_x*nbr_noeud_z
        DO j=1,nbr_noeud_x*nbr_noeud_z
            matrice_calcul_KM(i,j)=matrice_rigidite(i,j)+matrice_masse(i,j) ! Matrice K + M
        ENDDO
    ENDDO

        
    IF (choix_affichage == 1) THEN
        print*,"        t = 0 s"
        OPEN(unit=15,FILE="resultat.dat",ACTION="WRITE")
        DO w=1,nbr_noeud_x*nbr_noeud_z
            WRITE(15,'(*(g15.9,:,"' // ACHAR(9) // '"))')COORD(w,2),coord(w,3),matrice_U_N(w)
        ENDDO
        CALL system(" gnuplot -p data.plt ")
        CLOSE(15)
    ENDIF

    entree=0 !Compteur pour l'écriture dans EXCEL
    temps_x=temps_experience/pas_temps
    n=1
    CALL execute_command_line('open musique_dattente.mp3 ')
    OPEN(unit=15,FILE="resultat.dat",ACTION="WRITE")

    DO i=0,temps_x
        CALL execute_command_line('clear')
        PRINT*,"Calcul en cours ."
        !Création de la Matrice B
        matrice_b=0
        DO  balayage_triangle = 1,nombre_triangle
            k = connexe(balayage_triangle,2)
            l = connexe(balayage_triangle,3)
            m = connexe(balayage_triangle,4)
            x_k = coord(k,2)
            z_k = coord(k,3)
            x_l = coord(l,2)
            z_l = coord(l,3)
            x_m = coord(m,2)
            z_m = coord(m,3)

            det_JE = abs((x_l-x_k)*(z_m-z_k)-(x_m-x_k)*(z_l-z_k))

            flux_k=calcul_flux(matrice_U_N(k),teneur_eau_sat,temperature_humide,weq,&
            &temperature_air_ambiant,vitesse_air_ambiant,masse_molaire_vap,cte_gaz_parfait,hr_air_ambiant)
            flux_l=calcul_flux(matrice_U_N(l),teneur_eau_sat,temperature_humide,weq,&
            &temperature_air_ambiant,vitesse_air_ambiant,masse_molaire_vap,cte_gaz_parfait,hr_air_ambiant)
            flux_m=calcul_flux(matrice_U_N(m),teneur_eau_sat,temperature_humide,weq,&
            &temperature_air_ambiant,vitesse_air_ambiant,masse_molaire_vap,cte_gaz_parfait,hr_air_ambiant)

            IF ( ( (connexe(balayage_triangle,5)) ==  1 ) .AND. (connexe(balayage_triangle,6) == 1)  ) THEN ! Si SOMMET K est sur une bordure et  Si SOMMET L est sur une bordure
                    !TRIANGLE AVEC KL en bordure 
                    matrice_b(k)=matrice_b(k)+ (-pas_temps * (sqrt( (x_l-x_k)**2 + (z_l-z_k)**2 )) *&
                    & ( (flux_k/(3*masse_volumique))+ (flux_l/(6*masse_volumique)) ) ) 
                    matrice_b(l)=matrice_b(l) + (-pas_temps * (sqrt( (x_l-x_k)**2 + (z_l-z_k)**2 ))*&
                    & ( (flux_k/(6*masse_volumique))+ (flux_l/(3*masse_volumique)) ) )
            ENDIF
            IF ( ( (connexe(balayage_triangle,5)) ==  1 ) .AND. (connexe(balayage_triangle,7) == 1)  ) THEN! Si SOMMET K est sur une bordure et Si SOMMET M est sur une bordure
                    !TRIANGLE AVEC KM en bordure
                    matrice_b(k)=matrice_b(k)+ (-pas_temps * (sqrt( (x_m-x_k)**2 + (z_m-z_k)**2 )) * &
                    & ( (flux_k/(3*masse_volumique))+ (flux_l/(6*masse_volumique))) ) 
                    matrice_b(m)=matrice_b(m) + (-pas_temps * (sqrt( (x_m-x_k)**2 + (z_m-z_k)**2 ))* &
                    & ( (flux_k/(6*masse_volumique))+ (flux_m/(3*masse_volumique))) )
            ENDIF
            IF ( ((connexe(balayage_triangle,6)) ==  1) .AND. (connexe(balayage_triangle,7) == 1 ) ) THEN ! Si SOMMET L est sur une bordure et Si SOMMET M est sur une bordure
                    !TRIANGLE AVEC LM en bordure
                    matrice_b(l)=matrice_b(l)+ (-pas_temps * (sqrt( (x_l-x_m)**2 + (z_l-z_m)**2 )) *&
                    & ((flux_l/(3*masse_volumique))+ (flux_m/(6*masse_volumique)) ) )
                    matrice_b(m)=matrice_b(m) + (-pas_temps * (sqrt( (x_l-x_m)**2 + (z_l-z_m)**2  ))*&
                    & ( (flux_l/(6*masse_volumique))+ (flux_m/(3*masse_volumique)) ))
            ENDIF
        ENDDO
        
        CALL execute_command_line('clear')
        PRINT*,"Calcul en cours .."
        OPEN(unit=125,FILE="matrice_b.txt",ACTION="WRITE")
        DO w=1,nbr_noeud_x*nbr_noeud_z
            WRITE(125,'(*(g15.9,:,"' // ACHAR(9) // '"))')coord(w,:),matrice_b(w)
        ENDDO
        !Ecriture tout les temps d'écriture
        IF(i == n*temps_ecriture)THEN
            n=n+1
            CALL ecriture(nbr_noeud_x*nbr_noeud_z,nbr_noeud_x,nbr_noeud_z,i,matrice_U_N,entree)
            IF (choix_affichage == 1) THEN
                print*,"EVOLUTION EN FONCTION DU TEMPS  :"
                print*,"       t =",i,"s"
                OPEN(unit=15,FILE="resultat.dat",ACTION="WRITE")
                DO w=1,nbr_noeud_x*nbr_noeud_z
                    WRITE(15,'(*(g15.9,:,"' // ACHAR(9) // '"))')COORD(w,2),coord(w,3),matrice_U_N(w)
                ENDDO
                CALL system(" gnuplot -p data.plt ") 
                CLOSE(15)
            ENDIF
        ENDIF

        second_membre=matmul(matrice_masse,matrice_U_N)+matrice_b  ! MUn + B

        IF (choix_methode == 1) THEN !METHODE DE GAUSS
            !Creation de Ag pour resolution de Gauss
            DO w =1,nbr_noeud_x*nbr_noeud_z
                Ag(w,nbr_noeud_x*nbr_noeud_z+1)=second_membre(w)
                DO j=1,nbr_noeud_x*nbr_noeud_z
                    Ag(w,j)=matrice_calcul_KM(w,j)
                ENDDO
            ENDDO
    
            CALL calcul_gauss(nbr_noeud_z*nbr_noeud_x,ag,matrice_U_N1)

        ELSEIF (choix_methode ==2)THEN ! DÉCOMPOSITION LU
            call solve_linear_system_lu(nbr_noeud_z*nbr_noeud_x,matrice_calcul_KM,second_membre,matrice_U_N1)
        ELSEIF (choix_methode == 3)THEN !CHOLESKY
            call cholesky_solve(matrice_calcul_KM,second_membre,matrice_U_N1,nbr_noeud_x*nbr_noeud_z)
        ENDIF
        
        matrice_U_N=(matrice_U_N1) 
        
        CALL execute_command_line('clear')
        PRINT*,"Calcul en cours ..."
    ENDDO

    CALL execute_command_line('clear')

    PRINT*," VOULEZ VOUS AFFICHEZ LA DERNIÈRE ITÉRATION ?"
    PRINT*,"        OUI ========>>>>>> 1"
    PRINT*,"        NON ========>>>>>> 2"
    READ(*,*)choix_affichage_final

    IF (choix_affichage_final == 1) THEN
        OPEN(unit=15,FILE="resultat.dat",ACTION="WRITE")
        DO w=1,nbr_noeud_x*nbr_noeud_z
            WRITE(15,'(*(g15.9,:,"' // ACHAR(9) // '"))')COORD(w,2),coord(w,3),matrice_U_N(w)
        ENDDO
        CLOSE(15)
        CALL system(" gnuplot -p data.plt ")
    ENDIF

    !CALL execute_command_line(" open resultat.xls ")
    !CALL execute_command_line(" open resultat.txt ")
    !CALL execute_command_line(" open resultat_COUTURE.txt ")

    CALL execute_command_line(" open photo_final.jpg")



    PRINT*,"RESULTATS DANS LES FICHIERS : "
    PRINT*,"    - resultat.xls"
    PRINT*,"    - resultat.txt"
    PRINT*,"    - resultat_COUTURE.txt"
    PRINT*,"    - resultat.dat"
    PRINT*,""

END PROGRAM cas2

SUBROUTINE lecture1(epaisseur,alpha, teneur_eau_ini, teneur_eau_sat, masse_volumique)
    DOUBLE PRECISION,INTENT(OUT) :: epaisseur, alpha, teneur_eau_ini, teneur_eau_sat, masse_volumique
    OPEN(unit=10,FILE="donnee.txt",ACTION="READ")
    READ(10,*)
    READ(10,*)epaisseur
    READ(10,*)alpha
    READ(10,*)teneur_eau_ini
    READ(10,*)teneur_eau_sat
    READ(10,*)masse_volumique
    CLOSE(10)
END SUBROUTINE lecture1

SUBROUTINE lecture2(vecteur_gravite,masse_molaire_vap, cte_gaz_parfait, vitesse_air_ambiant)
    DOUBLE PRECISION,INTENT(OUT) :: vecteur_gravite, masse_molaire_vap, cte_gaz_parfait, vitesse_air_ambiant
    DOUBLE PRECISION :: bin
    INTEGER :: i
        OPEN(unit=10,FILE="donnee.txt",ACTION="READ")
        READ(10,*)
        DO i=1,5
            READ(10,*)bin
        ENDDO
        READ(10,*)vecteur_gravite
        READ(10,*)masse_molaire_vap
        READ(10,*)cte_gaz_parfait
        READ(10,*)vitesse_air_ambiant
        CLOSE(10)
END SUBROUTINE lecture2

SUBROUTINE lecture3(temperature_air_ambiant, hr_air_ambiant,Dzz,pas_temps,temps_experience,beta,nbr_noeud_x,&
    & nbr_noeud_z,temperature_humide,temps_ecriture)
    DOUBLE PRECISION,INTENT(OUT) ::temperature_air_ambiant, hr_air_ambiant,Dzz,pas_temps ,beta,temperature_humide
    INTEGER, INTENT(OUT):: nbr_noeud_x,nbr_noeud_z,temps_experience,temps_ecriture
    DOUBLE PRECISION :: bin
        OPEN(unit=10,FILE="donnee.txt",ACTION="READ")
        READ(10,*)
        DO i=1,9
            READ(10,*)bin
        ENDDO
        READ(10,*)temperature_air_ambiant
        READ(10,*)hr_air_ambiant
        READ(10,*)dzz
        READ(10,*)pas_temps
        READ(10,*)temps_experience
        READ(10,*)temps_ecriture
        READ(10,*)beta
        READ(10,*)nbr_noeud_x
        READ(10,*)nbr_noeud_z
        READ(10,*)temperature_humide
        CLOSE(10)
endsubroutine lecture3

function calcul_flux(teneur_eau,teneur_eau_sat,temperature_humide,weq,temperature_air_ambiant,&
    &vitesse_air_ambiant,masse_molaire_vap,cte_gaz_parfait,hr_air_ambiant) result(flux)
    DOUBLE PRECISION,INTENT(in) :: teneur_eau,teneur_eau_sat,temperature_humide,temperature_air_ambiant,weq
    DOUBLE PRECISION, INTENT(in) :: vitesse_air_ambiant,masse_molaire_vap,cte_gaz_parfait,hr_air_ambiant
    DOUBLE PRECISION :: pvsat,temperature,rhov,calcul_pvsat
    DOUBLE PRECISION :: a,b,at,bt,ct,aw,hm,pv_sat_air_ambiant,rho_inf
    DOUBLE PRECISION :: flux
    
    IF (teneur_eau > teneur_eau_sat) THEN
        temperature=temperature_humide+273.15
        pvsat=calcul_pvsat(temperature-273.15)
        rhov=(masse_molaire_vap*pvsat)/(cte_gaz_parfait*(temperature))
    ELSE
        a=(temperature_air_ambiant-temperature_humide)/(temperature_air_ambiant*temperature_humide*(teneur_eau_sat-weq))
        b= (1/temperature_humide) - a*teneur_eau_sat
        temperature=1/(a*teneur_eau+b) + 273.15 ! MISE EN K
        at = (-2.86d-5)*(temperature**2)-(1.07d-2)*temperature+10.24
        bt = (-5.41d-4)*temperature+1.01
        ct = (4.97d-6)*temperature**2-(2.67d-3)*temperature+0.35
        aw = exp(-at*((bt)**(100*teneur_eau))+ct)
        pvsat=calcul_pvsat(temperature-273.15)
        rhov=(masse_molaire_vap * (aw*pvsat) )/(cte_gaz_parfait*(temperature))
    END IF
    
    hm = (9.454d-3)*(vitesse_air_ambiant**(0.5003)) 

    pv_sat_air_ambiant = calcul_pvsat(temperature_air_ambiant)

    rho_inf = (masse_molaire_vap* hr_air_ambiant * pv_sat_air_ambiant )/(cte_gaz_parfait*(273.15+temperature_air_ambiant))

    flux=hm*(rhov-rho_inf)

END function calcul_flux

function calcul_pvsat(temperature) result(pvsat_result)
    DOUBLE PRECISION :: temperature ! en °C
    DOUBLE PRECISION :: pvsat_result

    pvsat_result=100.*exp(75.051-7293.3/(273.15+temperature)-8.593*&
    &log((273.15+temperature))+0.00617*(273.15+temperature))
    
end function calcul_pvsat

SUBROUTINE ecriture(n,nx,nz,temps,w,entree) 
    IMPLICIT NONE
    DOUBLE PRECISION, DIMENSION(n) :: w
    INTEGER :: i,temps,n,entree,nz,nx,j

    OPEN(unit=10,FILE="resultat.xls",ACTION="WRITE")
    entree=entree+1
    IF (entree == 1) THEN
        WRITE(10,'(*(g15.9,:,"' // ACHAR(9) // '"))')'temperature (en s)','Humidite'
        WRITE(10,'(*(g15.9,:,"' // ACHAR(9) // '"))')temps,(W(j),j=i,(nx-1)*nz+i,nz)
    ENDIF
    DO i=1,nz
        WRITE(10,'(*(g15.9,:,"' // ACHAR(9) // '"))')temps,(W(j),j=i,(nx-1)*nz+i,nz)
    ENDDO
    WRITE(10,'(*(g15.9,:,"' // ACHAR(9) // '"))')


    OPEN(unit=11,FILE="resultat.txt",ACTION="WRITE")
    OPEN(unit=12,FILE="resultat_COUTURE.txt",ACTION="WRITE")
    Do i=1,nz
        write(11,'(*(g15.9,:,"' // ACHAR(9) // '"))') (W(j),j=i,(nx-1)*nz+i,nz)
        write(12,'(500(a,g15.9))') (',',W(j),j=i,(nx-1)*nz+i,nz)
    enddo
    WRITE(11,'(*(g15.9,:,"' // ACHAR(9) // '"))')
    WRITE(12,'(*(g15.9,:,"' // ACHAR(9) // '"))')
    

END SUBROUTINE ecriture

FUNCTION calcul_th(Th,t_ambiant,hr_ambiant) result(f_sortie)
    DOUBLE PRECISION :: th,calcul_pvsat,pvsat_th,pvsat_t,f_sortie,t_ambiant,hr_ambiant

    pvsat_th = calcul_pvsat(Th)
    pvsat_t = calcul_pvsat(t_ambiant)
    !FONCTION : (hm/ht).(Mv/R).[ (Pv(th)/Th) - HR*(Pvsat(T)/T) ]
    f_sortie = (1.0/1000)*(0.018/8.3143)*( (pvsat_th)/(Th+273) - ((hr_ambiant*pvsat_t)/(t_ambiant+273.15)) )* &
    &(2185*(Th)+2501d3) + (Th - t_ambiant)

END FUNCTION calcul_th

SUBROUTINE dichotomie_Th(f,x,a,b,hr_ambiant,t_ambiant)
    IMPLICIT NONE

    DOUBLE PRECISION, EXTERNAL :: f
    DOUBLE PRECISION, INTENT(INOUT) :: hr_ambiant,t_ambiant
    DOUBLE PRECISION, INTENT(OUT):: x
    DOUBLE PRECISION :: a,b,eps

    eps=1e-8
    DO WHILE (abs(b-a)>eps)
        x=(b+a)/2
        IF (f(a,t_ambiant,hr_ambiant)*f(x,t_ambiant,hr_ambiant) > 0) THEN
            a=x
        ELSE
            b=x
        ENDIF
    END DO
END SUBROUTINE dichotomie_Th

!https://quasiengineer.dev/tech/engg/gauss-elimination-in-fortran/
SUBROUTINE calcul_gauss(n,ag,x)
    use, intrinsic :: iso_fortran_env, only : dp=>real64, error_unit, input_unit
    INTEGER :: n, pos ! n - Number of unknowns
    real(dp) :: ag(n,n+1), x(n)! Augumented matrix and solution vector
    real(dp) :: factor
    
    !Forward elimination
    DO k = 1, n-1
        ! Partial scaled pivoting
        pos = maxloc( abs(ag(k:n, k)/maxval(ag(k:n, :), dim=2)), dim=1 )
        j = k + pos - 1
        IF (.not. j == k) CALL swap(ag(j,:), ag(k,:))
        ! Elimination
        DO i = k+1, n
            factor = ag(i,k) / ag(k,k)
            ag(i, k:) = ag(i, k:) - factor*ag(k, k:)
        END DO
    END DO

    !Back substitution
    x(n) = ag(n, n+1) / ag(n,n)
    DO i = n-1, 1, -1
        x(i) = ( ag(i,n+1) - dot_product(ag(i, i+1:n), x(i+1:n)) ) / ag(i,i)
    END DO

    contains
    elemental SUBROUTINE swap(a, b)
        DOUBLE PRECISION, INTENT(INOUT) :: a, b
        DOUBLE PRECISION :: temp
        temp = a
        a = b
        b = temp
    END SUBROUTINE swap
    
END SUBROUTINE calcul_gauss

!CAS 1 - Stéphane GIBOUT
!MODELISATION THERMODYNAMIQUE I - Elodie LE GUEN & Erwan BENHENOU
SUBROUTINE DECOMPO_LU(matrice_entree,n,matrice_l,matrice_u)
    INTEGER,INTENT(IN):: n
    DOUBLE PRECISION, DIMENSION(n,n),INTENT(IN) ::matrice_entree
    DOUBLE PRECISION, DIMENSION(n,n),INTENT(OUT) :: matrice_l,matrice_u

    INTEGER :: i,k,j
    DOUBLE PRECISION :: somme_u,somme_l

    DO i =1,n
        DO j=1,n
            matrice_u(i,j)=0
            IF (i == j) THEN
                matrice_l(i,j)=1
            ELSE
                matrice_l(i,j)=0
            ENDIF
        ENDDO
    ENDDO

    DO j=1,n
        somme_u=0
        somme_l=0 
        DO i=1,n
            IF (i > j) THEN
                DO k=1,j-1
                    somme_l = somme_l+ matrice_l(i,k)* matrice_u(k,j)
                ENDDO
                matrice_l(i,j)=(1/(matrice_u(j,j))) * (matrice_entree(i,j) - somme_l )
            ELSE IF (i <= j) THEN
                DO k=1,i-1
                    somme_u = somme_u+ matrice_l(i,k)*matrice_u(k,j)
                ENDDO
                matrice_u(i,j)=matrice_entree(i,j) - somme_u
            ENDIF
        ENDDO
    ENDDO
    
END SUBROUTINE DECOMPO_LU

!CHATGPT
subroutine solve_linear_system_lu(n, A, B, X)
  
    implicit none
    
    integer, intent(in) :: n
    real*8, intent(inout) :: A(n,n)
    real*8, intent(inout) :: B(n)
    real*8, intent(out) :: X(n)
    real*8 :: L(n,n), U(n,n)
    integer :: i, j, k
    
    ! LU factorization using Gaussian elimination
    do k = 1, n-1
        do i = k+1, n
            L(i,k) = A(i,k)/A(k,k)
            do j = k+1, n
                A(i,j) = A(i,j) - L(i,k)*A(k,j)
            end do
            A(i,k) = L(i,k)
        end do
    end do
    U = A
    
    ! forward substitution to solve LY=B for Y
    do i = 1, n
        X(i) = B(i)
        do j = 1, i-1
            X(i) = X(i) - L(i,j)*X(j)
        end do
    end do
    
    ! backward substitution to solve UX=Y for X
    do i = n, 1, -1
        do j = i+1, n
            X(i) = X(i) - U(i,j)*X(j)
        end do
        X(i) = X(i)/U(i,i)
    end do
  
end subroutine solve_linear_system_lu

subroutine cholesky_solve(A, B, X, N)
    implicit none
    integer, intent(in) :: N
    real(kind=8), intent(in) :: A(N,N), B(N)
    real(kind=8), intent(out) :: X(N)
    real(kind=8) :: L(N,N), Y(N), Z(N) ,one(n,n)
    integer :: i, j, k
  
    DO i =1,N
        Do j=1,N
            one(i,j)=1.0
        enddo
    enddo
    ! Compute the Cholesky decomposition of A
    L(1,1) = sqrt(A(1,1))
    do j = 2, N
      L(j,1) = A(j,1) / L(1,1)
    end do
    do i = 2, N
      L(i,i) = sqrt(A(i,i) - dot_product(L(i,1:i-1)**2, one(i,1:i-1)) )
      do j = i+1, N
        L(j,i) = (A(j,i) - dot_product(L(j,1:i-1)*L(i,1:i-1), one(i,1:i-1) ) ) / L(i,i)
      end do
    end do
  
    ! Solve the triangular system L*Y=B using forward substitution
    Y(1) = B(1) / L(1,1)
    do i = 2, N
      Y(i) = (B(i) - dot_product(L(i,1:i-1), Y(1:i-1))) / L(i,i)
    end do
  
    ! Solve the triangular system L^T*X=Y using backward substitution
    Z = 0.0
    do i = 1, N
      X(i) = Y(i) / L(i,i)
      do j = i+1, N
        Z(j) = Z(j) + L(j,i)*X(i)
      end do
    end do
    do i = N, 1, -1
      X(i) = (X(i) - dot_product(L(i,i+1:N), Z(i+1:N))) / L(i,i)
    end do
  
end subroutine cholesky_solve  