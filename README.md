# Cinetica_Quimica

> Status: Projeto em Desenvolvimento

Este projeto está sendo densenvolvido por [Ana Clara Loponi](https://github.com/AnaLoponi) e [João Pedro Aroucha](https://github.com/jpab2004). O repositório contém o código, material textual produzido e informações relacionadas ao projeto de simulação de reações químicas, desenvolvido como parte da disciplina ___Cinética Química___, orientada pelo Profº [Amauri Jardim de Paula](http://lattes.cnpq.br/3912955675367941), na ___Ilum - Escola de Ciência___.

## Tópicos Centrais do Problema

A cinética química é a área da química que trata das taxas (velocidades) nas quais as reações químicas ocorrem e dos fatores que influenciam essas taxas. Este projeto visa simular uma reação química utilizando os conceitos e equações da cinética química. Os tópicos centrais abordados no projeto incluem:

1. **Mecanismos de Reação:** Exploração das etapas individuais e espécies intermediárias em reações químicas, especialmente as complexas que ocorrem através de etapas elementares que formam uma reação global.

2. **Leis de Velocidade:** Análise da relação entre a taxa de reação e as concentrações dos reagentes, representadas por leis matemáticas que expressam essa relação.

3. **Energia de Ativação:** Compreensão da barreira de energia que os reagentes devem superar para ocorrer a conversão entre reagentes e produtos, e a influência da energia de ativação na dependência da temperatura da reação.

4. **Ordem da Reação:** Investigação de como a concentração de cada reagente afeta a taxa de reação e como a ordem da reação determina essa dependência.

5. **Fatores que Afetam a Taxa de Reação:** Exploração de fatores externos como temperatura, pressão e catalisadores, que podem alterar significativamente a velocidade de uma reação química.

## Objetivos

O objetivo principal deste projeto é criar um modelo computacional de simulação de cinética química que leve em consideração os diversos conceitos e parâmetros relacionados à reação química. Os objetivos específicos incluem:

1. Definir as condições iniciais do sistema com base na reação química escolhida.

2. Estabelecer as velocidades das partículas reagentes utilizando uma função estatística apropriada.

3. Identificar as leis de velocidade específicas para a reação química selecionada.

4. Considerar as condições de energia de ativação relevantes para a reação química.

5. Desenvolver um modelo computacional em Python que simule a cinética da reação química.

## Desafios

Durante o semestre, este projeto enfrentará três desafios significativos que serão abordados em etapas sucessivas, visando à construção do modelo de simulação. 

### Desafio 1: Descrição do Movimento de Partículas
> <sub> Em desenvolvimento. [Desafio 1](https://github.com/jpab2004/Cinetica_Quimica/blob/main/main.py)</sub>

Nesta etapa, desenvolvemos um modelo inicial que descreve as variáveis físicas relevantes para cada partícula (átomo ou molécula) envolvida na reação química. A descrição inclui informações de posição e velocidade das partículas em 2 e 3 dimensões. Este modelo está obedecendo aos princípios de conservação do momento e da energia, considerando também a natureza estatística das partículas e sua dependência de variáveis termodinâmicas, como a temperatura. As partículas estão sendo consideradas como gases ideais.

### Desafio 2: Definição da Reação Química e Leis de Velocidade
> <sub> Não implementado. </sub>

Nesta fase, iremos incorporar a conversão dos reagentes em produtos, caracterizando a reação química. Serão consideradas as leis de velocidade de conversão, equações diferenciais que dependem das concentrações dos reagentes e da ordem da reação. A experimentação será essencial para definir essas leis de velocidade.

### Desafio 3: Estereoquímica e Energia de Ativação
> <sub> Não implementado. </sub>

No último desafio, será considerado o papel da estereoquímica na formação dos produtos, levando em conta os aspectos espaciais e energéticos da conversão dos reagentes. Também será abordado o efeito dos catalisadores na alteração das leis de velocidade. Questões estereoquímicas, energéticas e catalíticas serão incorporadas ao modelo.

## Estrutura do Repositório - Guia

* main.py
> <sub> Código implementado. [main.py](https://github.com/jpab2004/Cinetica_Quimica/blob/main/main.py)</sub>
* LICENSE
> <sub>Necessário para licenciar o uso do repositório. Versão: GNU General Public License v3.0. [Saiba Mais](https://docs.github.com/en/communities/setting-up-your-project-for-healthy-contributions/adding-a-license-to-a-repository)</sub>
* .gitignore
> <sub> Arquivo texto que indica regras ao Git sobre quais arquivos ou pastas ele deve ignorar em um projeto. [Saiba Mais](https://docs.github.com/en/get-started/getting-started-with-git/ignoring-files)</sub>
* README.md
> <sub> Onde estão armazenadas todas as informações que você acabou de ler ✔️</sub>

**Importante:** O conteúdo deste repositório é voltado para fins educacionais e de aprendizado no âmbito do projeto de simulação de cinética química.
