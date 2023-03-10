# fisica-do-clima

Repositório criado para guardar os scripts utilizados na confecção do seminário da disciplina Física do Clima.
Os dados utilizados não estão inclusos

## Objetivo

 O objetivo deste trabalho é realizar uma análise climatológica das ciclogêneses no Oceano Atlântico Sudoeste (OAS) e estudar a sua relação com a TSM e fluxos de calor latente e calor sensível, utilizando como fonte dados da reanálise ERA-5 para um novo período climatológico entre 1991 a 2021, até então ainda não explorado.

## Metodologia
 
 As ciclogênese serão identificadas a partir de um algoritmo de tracking que utiliza como input dados de vorticidade relativa em baixos níveis do ERA5 (6 em 6 horas). Serão considerados ciclones todo sistema que possui duração igual ou superior a 48 horas (Hodges et al., 2002).
 De acordo com os resultados da analise de sensibilidade para o ano de 2021 (validado com cartas sinoticas), sistemas que apresentarem vorticidade relativa superior a -4 x 10^(-5) 1/s após 60 horas da sua identificação pelo algoritmo também serão filtrados. 