import telebot
import config
import contents
#import dbworker
import BioMethods
import os
import requests


bot = telebot.TeleBot(config.TOKEN, parse_mode='HTML')

#СТАРТОВЫЕ СООБЩЕНИЯ
@bot.message_handler(commands=['start'])
def start_message(message):
    bot.send_message(message.chat.id, contents.base_menu_message, reply_markup=contents.menu_base)
    bot.send_message(message.chat.id, contents.design_menu_message, reply_markup=contents.menu_design)
    bot.send_message(message.chat.id, contents.file_menu_message, reply_markup=contents.menu_file)

#Callback Querry









if __name__ == "__main__":
    bot.polling()
