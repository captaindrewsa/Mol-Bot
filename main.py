import telebot
import config
#import contents
#import dbworker
import BioMethods
import os
import requests


bot = telebot.TeleBot(config.TOKEN, parse_mode='HTML')










if __name__ == "__main__":
    bot.polling()
