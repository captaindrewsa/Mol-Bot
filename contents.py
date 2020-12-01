import re
import telebot

#Кнопочки
BTNS = [
    telebot.types.InlineKeyboardButton("Разворот", callback_data='reverse'),
    telebot.types.InlineKeyboardButton("Комплиментарность", callback_data='compl'),
    telebot.types.InlineKeyboardButton("Температура плавления олигонуклеотида", callback_data='temp_olig'),
    telebot.types.InlineKeyboardButton("Lambda-red Recombination", callback_data='lrr'),
    telebot.types.InlineKeyboardButton("Создание .fasta файла", callback_data='create_fasta'),
    telebot.types.InlineKeyboardButton("БОЛЬШОЙ ТЕКСТ МНОГО БУКАВ", callback_data='tuc'),
    telebot.types.InlineKeyboardButton("БОЛЬШОЙ ТЕКСТ МНОГО БУКАВ", callback_data='tuc')]

#Месседжи
base_menu_message = '''Тут будут <b>базовые</b> функции'''
design_menu_message = '''Тут будет меню <b>экспериментов</b>'''
file_menu_message = '''Тут будет функционал с <b>файлами</b>'''


help_message = '''На данный момент бот умеет следующее:\n
•Работать с файлами формата <b>fasta</b>
•Разворачивать нуклеотидную последовательность
•Считать температуру плавления олигов
•Выдавать комплементарную цепочку (RNA и DNA)\n  
Чтобы вызвать меню воможностей - введите /menu
'''
#Менюшки
menu_base = telebot.types.InlineKeyboardMarkup()
menu_base.row(BTNS[0]).row(BTNS[1]).row(BTNS[2])

menu_design = telebot.types.InlineKeyboardMarkup()
menu_design.row(BTNS[3])

menu_file = telebot.types.InlineKeyboardMarkup()
menu_file.row(BTNS[4])

menu_all = telebot.types.ReplyKeyboardMarkup(resize_keyboard=True, one_time_keyboard=True)
menu_all.row('Разворот','Комплементарность', 'Температура плавления олига') #Меню полное

menu_compl = telebot.types.ReplyKeyboardMarkup(resize_keyboard=True, one_time_keyboard=True)
menu_compl.row('Комплементарная DNA','Комплементарная RNA') #Меню выбора комплиментарности

menu_menu = telebot.types.ReplyKeyboardMarkup(resize_keyboard=True, one_time_keyboard=True)
menu_menu.row('/menu','/fasta') #Меню после функции

menu_tempolig = telebot.types.ReplyKeyboardMarkup(resize_keyboard=True, one_time_keyboard=True)
menu_tempolig.row('DNA','RNA') #Меню выбора олига
'''
menu_file = telebot.types.ReplyKeyboardMarkup(resize_keyboard=True, one_time_keyboard=True)
menu_file.row('Пункт 1','Пункт 2', "Пункт...") #Меню выбора олига
'''
menu_all_na = telebot.types.ReplyKeyboardMarkup(resize_keyboard=True, one_time_keyboard=True)
menu_all_na.row('/все')










